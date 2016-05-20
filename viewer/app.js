var express = require('express');
var app = express();
var server = require('http').Server(app);
var io = require('socket.io')(server);
var browserify = require('browserify-middleware');

var declared = function(v) {
	return (typeof v != 'undefined');
};

app.get('/', function(req, res) {
	res.sendFile(__dirname + '/public/index.html');
});
app.get('/js/main.js', browserify(__dirname + '/public/js/main.js'));

app.get(/^(.+)$/, function(req, res) {
	res.sendFile(__dirname + '/public/' + req.params[0]);
});

var port = process.argv.indexOf("-p") != -1 ? process.argv[process.argv.indexOf("-p") + 1] : 3000;
server.listen(port, '0.0.0.0', function() {
	console.log('Grgen viewer running on port ' + port);
});


process.stdin.setEncoding('utf8');

var GRNFrames = [];

io.on('connection', function(socket) {
	console.log("new socket connection");
	GRNFrames.forEach(function(f) {
		console.log("sending frame");
		socket.emit('GRNFrame', f);
	});
});

var bracketsLvl = 0;
process.stdin.on('readable', function() {
	var chunk = process.stdin.read();
	if (chunk !== null) {
		var buf = "";
		while (chunk.length > 0) {
			var nextChar = chunk[0];
			chunk = chunk.substring(1);
			if (nextChar == "{")
				bracketsLvl++;
			else
			if (nextChar == "}")
				bracketsLvl--;
			buf += nextChar;
			if (bracketsLvl === 0) {
				if (buf.length > 1) {
					try {
						data = JSON.parse(buf);
						if (declared(data.GRN)) {
							data.time = new Date().getTime();
							data.frameId = GRNFrames.length;
							GRNFrames.push(data.GRN);
							io.sockets.emit('GRNFrame', data.GRN);
						}
					} catch (e) {
						console.log('ERROR:' + e);
						console.log('chunk : ' + buf);
					}
				}
				buf = '';
			}
		}
	}
});

process.stdin.on('end', function() {
	process.stdout.write('end');
});
