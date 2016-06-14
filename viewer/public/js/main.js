var tools = require('./tools.js');
tools.addEqualsToArray();
var Victor = require('victor');
var io = require('socket.io-client');
var Chart = require('chart.js');
var $ = require('jquery-browserify');
var res = 2;
var rdmColor = require('randomColor');


function range(start, end) {
	var foo = [];
	for (var i = start; i <= end; i++) {
		foo.push(i);
	}
	return foo;
}

var plotLines = {};
var data = {
	labels: []
};

function getDataRange() {
	var max = 0;
	for (var key in plotLines)
		if (plotLines.hasOwnProperty(key) && plotLines[key].values.length > max) max = plotLines[key].values.length;
	return range(0, max - 1);
}


var options = {
	maintainAspectRatio: false,
	animation: false,
	scales: {
		xAxes: [{
			ticks: {
				display: false
			},
			gridLines: {
				display: false
			}
		}]
	},
	hover: {
		mode: 'label'
	},
	width: 900,
	responsive: false
};


PIXI.Point = Victor;

var renderer = PIXI.autoDetectRenderer(window.innerWidth, window.innerHeight, {
	antialias: true,
	resolution: res,
	autoResize: true
});


var stage = new PIXI.Container();
stage.interactive = true;
stage.hitArea = new PIXI.Rectangle(0, 0, window.innerWidth, window.innerHeight);
renderer.view.style.position = "absolute";
renderer.view.style.display = "block";

function resize() {
	renderer.resize(window.innerWidth, window.innerHeight);
}

document.body.appendChild(renderer.view);
window.addEventListener("resize", resize);
var concentrationPlotCtx = document.getElementById("concentration_plot").getContext("2d");
concentrationPlotCtx.canvas.width = $("#concentration_plot").parent().width();
var concentrationChart = new Chart(concentrationPlotCtx, {
	type: 'line',
	data: data,
	scales: {
		xAxes: [{
			ticks: {
				autoSkip: true
			}
		}]
	},
	options: options
});

function updatePlotData() {
	data.labels = getDataRange();
	data.datasets = [];
	for (var key in plotLines) {
		if (plotLines.hasOwnProperty(key)) {
			var line = {
				label: key,
				fill: false,
				borderWidth: 1,
				borderColor: "#" + plotLines[key].color.toString(16),
				backgroundColor: "#" + plotLines[key].color.toString(16),
				lineTension: 0.0,
				pointRadius: 0,
				data: plotLines[key].values
			};
			data.datasets.push(line);
		}
	}
	concentrationChart.update();
}



var PROT_SIZE = 30;
var PROT_INNER_SIZE = 10;
var PROT_BASE_POSITION = new Victor(window.innerWidth / (2 * res), window.innerHeight / (2 * res));
var PROT_MIN_DIST = 150;
var PROT_MAX_DIST = 400;
var PROT_INPUT = new Victor(50, window.innerHeight / (2 * res));
var PROT_OUTPUT = new Victor(window.innerWidth / res - 50, window.innerHeight / (2 * res));
var SPRINGS_K = 1;
var SPRINGS_C = 0.15;
var VISCO = 0.5;

var grnscene = new PIXI.Container();
grnscene.hitArea = new PIXI.Rectangle(0, 0, window.innerWidth, window.innerHeight);

function Protein(opt) {
	// PROTEIN CLASS //

	PIXI.Container.apply(this);

	// ---------------------
	//         Guts
	// ---------------------
	this.c = tools.declared(opt) && tools.declared(opt.c) ? opt.c : 0;
	this.concentrationHistory = [this.c];
	this.prevc = tools.declared(opt) && tools.declared(opt.prevc) ? opt.prevc : 0;
	this.coords = tools.declared(opt) && tools.declared(opt.coords) ? opt.coords : [];
	this.name = tools.declared(opt) && tools.declared(opt.name) ? opt.name : '';
	this.input = tools.declared(opt) && tools.declared(opt.input) ? opt.input : false;
	this.output = tools.declared(opt) && tools.declared(opt.output) ? opt.output : false;

	// ---------------------
	//         Skin 
	// ---------------------
	this.position = new Victor(PROT_BASE_POSITION.x, PROT_BASE_POSITION.y);
	this.velocity = new Victor();
	this.force = new Victor();
	this.fixed = false;
	this.selected = false;
	this.dragged = false;

	this.color = tools.declared(opt) && tools.declared(opt.color) ? opt.color : parseInt(rdmColor({
		luminosity: 'light'
	}).substr(1), 16);
	this.circles = new PIXI.Graphics();
	this.text = new PIXI.Text(this.name, {
		font: '9px Helvetica',
		fill: 0xffffff,
		align: 'center'
	});
	this.text.resolution = 2;
	this.text.x = this.position.x - (this.text.width / 2);
	this.text.y = this.position.y - (this.text.height / 2);

	this.addChild(this.circles);
	this.addChild(this.text);

	// --------------------
	//        Update
	// --------------------
	this.updateConcentration = function(newC) {
		this.concentrationHistory.push(newC);
		this.c = newC;
	};

	this.update = function(dt) {
		if (!this.fixed && !this.dragged) {
			oldVel = this.velocity;
			this.force.subtract(oldVel.clone().multiplyScalar(VISCO));
			this.velocity.add(this.force.clone().multiplyScalar(dt));
			this.position.add(this.velocity.clone().add(oldVel).multiplyScalar(dt * 0.5));
		}
		this.force = new Victor();

		this.circles.clear();
		if (this.selected)
			this.circles.lineStyle(2, this.color); // (thickness, color)
		if (this.dragged)
			this.circles.lineStyle(2, 0xFFFFFF); // (thickness, color)
		if (!this.selected && !this.dragged)
			this.circles.lineStyle(0.7, 0x999999); // (thickness, color)
		this.circles.drawCircle(this.position.x, this.position.y, PROT_SIZE); // (x,y,radius)
		this.circles.lineStyle(0, 0x999999);
		this.circles.beginFill(this.color);
		this.circles.drawCircle(this.position.x, this.position.y, PROT_INNER_SIZE + (PROT_SIZE - PROT_INNER_SIZE) * this.c);
		this.circles.lineStyle(0.7, 0xBBBBBB);
		this.circles.beginFill(0x444444);
		this.circles.drawCircle(this.position.x, this.position.y, PROT_INNER_SIZE);
		this.circles.endFill();
		this.text.x = this.position.x - (this.text.width / 2);
		this.text.y = this.position.y - (this.text.height / 2);
	};
	this.interactive = true;
	this.hasBeenDragged = false;
	this.mousedown = function(event) {
		this.dragged = true;
		this.hasBeenDragged = false;
		event.stopPropagation();
	};
	this.mouseup = function(event) {
		this.dragged = false;
		if (!this.hasBeenDragged) {
			this.selected = !this.selected;
			if (this.selected) {
				plotLines[this.name] = {
					values: this.concentrationHistory,
					color: this.color
				};
				updatePlotData();
			} else {
				delete plotLines[this.name];
				updatePlotData();
			}
		}
		this.hasBeenDragged = false;
	};
	this.mousemove = function(event) {
		if (this.dragged) {
			this.hasBeenDragged = true;
			var p = event.data.getLocalPosition(grnscene);
			this.x = p.x / 2;
			this.y = p.y / 2;
		}
	};
}

Protein.prototype = new PIXI.Container();

function Spring(a, b, l, k, c) {
	this.A = a;
	this.B = b;
	this.K = k;
	this.C = c;
	this.L = l;
	this.currentL = l;
	this.update = function(dt) {
		prevL = this.currentL;
		dir = this.B.position.clone().subtract(this.A.position);
		this.currentL = dir.length();
		while (this.currentL === 0) {
			this.A.position.add(new Victor(Math.random() - 0.5, Math.random() - 0.5));
			dir = this.B.position.clone().subtract(this.A.position);
			this.currentL = dir.length();
		}
		dir.divideScalar(this.currentL);
		speed = (this.currentL - prevL) / dt;
		X = this.L - this.currentL;
		dir.multiplyScalar(this.K * X - speed * this.C);
		this.A.force.subtract(dir);
		this.B.force.add(dir);
	};
}

function GRN(opt) {
	this.params = tools.declared(opt) && tools.declared(opt.params) ? opt.params : [];
	this.proteins = tools.declared(opt) && tools.declared(opt.proteins) ? opt.proteins : [];
	this.subNets = tools.declared(opt) && tools.declared(opt.subNets) ? opt.subNets : [];
	this.springs = [];
	this.recreateSprings = function() {
		this.springs = [];
		for (var i = 0; i < this.proteins.length; ++i) {
			for (var j = i + 1; j < this.proteins.length; ++j) {
				this.springs.push(new Spring(this.proteins[i], this.proteins[j], tools.lerp(PROT_MIN_DIST, PROT_MAX_DIST, Math.abs(this.proteins[i].coords[0] - this.proteins[j].coords[0])), SPRINGS_K, SPRINGS_C));
			}
		}
		var nbI = 0;
		var nbO = 0;
		for (i = 0; i < this.proteins.length; ++i) {
			if (this.proteins[i].input && !this.proteins[i].output) {
				++nbI;
				this.proteins[i].fixed = true;
				this.proteins[i].position = new Victor(PROT_INPUT.x, PROT_INPUT.y + (nbI % 2 === 0 ? -1 : 1) * (Math.floor(nbI * 0.5) * 0.4 * PROT_MIN_DIST));
			} else if (!this.proteins[i].input && this.proteins[i].output) {
				++nbO;
				this.proteins[i].fixed = true;
				this.proteins[i].position = new Victor(PROT_OUTPUT.x, PROT_OUTPUT.y + (nbO % 2 === 0 ? -1 : 1) * (Math.floor(nbO * 0.5) * 0.4 * PROT_MIN_DIST));
			}
		}
	};
	this.update = function(dt) {
		for (var i = 0; i < this.springs.length; ++i) {
			this.springs[i].update(dt);
		}
		for (i = 0; i < this.proteins.length; ++i) {
			this.proteins[i].update(dt);
		}
	};
}

var indexToName = {};

function refreshIndexToName() {
	for (var index in fg.proteins.namedIn) {
		if (fg.proteins.namedIn.hasOwnProperty(index)) {
			indexToName[fg.proteins.namedIn[index]] = index;
		}
	}
	for (index in fg.proteins.namedOut) {
		if (fg.proteins.namedOut.hasOwnProperty(index)) {
			indexToName[fg.proteins.namedOut[index]] = index;
		}
	}
}
GRN.prototype.updateFromFrame = function(fg) {
	if (fg.proteins.plist.length !== this.proteins.length) {
		console.log("Creating protein array from scratch");
		this.proteins = fg.proteins.plist.map(function(p, i) {
			return new Protein({
				c: p.c,
				coords: p.coords,
				name: indexToName[i] !== undefined ? indexToName[i] : "Default " + i,
				input: p.I,
				output: p.O
			});
		});
		//TODO: remove old proteins from stage
		this.proteins.forEach(function(p) {
			grnscene.addChild(p);
		});
		this.recreateSprings();
		refreshIndexToName();
	} else {
		this.proteins.forEach(function(p, i) {
			if (!p.coords.equals(fg.proteins.plist[i].coords)) {
				console.log("different coords");
			}
			p.updateConcentration(fg.proteins.plist[i].c);
			name = indexToName[i] !== undefined ? indexToName[i] : "Default " + i;
		});
	}
	if (!this.params.equals(fg.params)) {
		console.log("different params");
		this.params = fg.params;
	}
	data.labels = getDataRange();
	concentrationChart.update();
};


var grn = new GRN();


var GRNFrames = [];
var socket = io.connect();
socket.on('connect', function(data) {
	console.log("connected to server");
	socket.on('GRNFrame', function(f) {
		console.log("received grn frame");
		if (!tools.declared(grn)) grn = new GRN();
		else {
			GRNFrames[f.frameID] = f;
			grn.updateFromFrame(f);
		}
	});
});


var FPS = new PIXI.Text('0 fps', {
	font: '9px Helvetica',
	fill: 0xaaaaaa,
});
FPS.resolution = 2;
FPS.x = 5;
FPS.y = 5;

stage.mousedown = stage.touchstart = function(data) {
	this.dragging = true;
};
stage.mouseup = stage.mouseupoutside = stage.touchend = stage.touchendoutside = function(data) {
	this.dragging = false;
};
stage.mousemove = stage.touchmove = function(data) {
	if (this.dragging) {
		var delta = new Victor(data.data.originalEvent.movementX, data.data.originalEvent.movementY);
		grnscene.x += delta.x;
		grnscene.y += delta.y;
	}
};

stage.addChild(FPS);
stage.addChild(grnscene);

var frameSinceLastFPSUpdt = 0;
var delta = 0;
lasttime = new Date().getTime();
requestAnimationFrame(animate);

function animate() {
	var currtime = new Date().getTime();
	var deltaSec = (currtime - lasttime) / 1000;
	delta += deltaSec;
	++frameSinceLastFPSUpdt;
	if (delta > 0.3) {
		FPS.text = Math.round(1 / (delta / frameSinceLastFPSUpdt)) + ' fps';
		frameSinceLastFPSUpdt = 0;
		delta = 0;
	}
	grn.update(Math.min(deltaSec, 0.1));
	renderer.render(stage);
	requestAnimationFrame(animate);
	lasttime = currtime;
}
