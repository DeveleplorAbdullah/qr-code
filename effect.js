// $(document).ready(function () {

// 	if(/Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent)){
// 		Webcam.set('constraints',{
// 			video: true,
// 			facingMode:  "environment"
// 		});
// 	}
	
// 		Webcam.attach('#example');
// 		setInterval(function(){ take_snapshot();}, 10);

// 	qrcode.callback = showInfo;

// });

function take_snapshot() {
	Webcam.snap(function (dataUrl) {
		qrCodeDecoder(dataUrl);
	});
}

// decode the img
function qrCodeDecoder(dataUrl) {
	qrcode.decode(dataUrl);
}

// show info from qr code
function showInfo(data) {
	$("#qrContent").text(data);
}

function callCamera(){
	console.log("call");
	if(/Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent)){
		Webcam.set('constraints',{
			video: true,
			facingMode:  "environment"
		});
	}
	
		Webcam.attach('#example');
		setInterval(function(){ take_snapshot();}, 10);

	qrcode.callback = showInfo;
}