var stepcounter = 0 ;
var maxsteps = 6;
var animate = false ;

var gapcost = 1 ;
var scorematch = 1;
var scoremismatch = -1;

function buildtable(){

	$('#dpmatrix').empty();
	$('#dpmatrix2').empty();

	// GCATGCU
	// GATTACA
	var a = $('input[name=seq1]').val();
	var b = $('input[name=seq2]').val();
	if (a.length > 0 && b.length >0){
		var alignment = Linear_alignment(a,b,1,-1,gapcost);
		console.log(alignment);
		for (i = 0; i <= b.length+1; i++) {
			$('#dpmatrix').append("<tr>");
			if(i<2){
				$('#dpmatrix').append("<th>-</th>");
			}else{
				$('#dpmatrix').append("<th>"+b.charAt(i-2)+"</th>");
			}
			for (j=0; j <= a.length; j++) {
				if(i===0){
					if(j===0){
						$('#dpmatrix').append("<th>-</th>");
					}else{
						$('#dpmatrix').append("<th>"+a.charAt(j-1)+"</th>");
					}
				}else{
					var num = alignment.Mat[j][i-1];
					var arrowid = "";
					if(alignment.Dir[j][i-1]===0){
						arrowid="diag";
					}
					else if(alignment.Dir[j][i-1]==-1){
						arrowid="left";
					}else if (alignment.Dir[j][i-1]==1) {
						arrowid="up";
					}else{
						arrowid = "noarrow";
					}

					var cell = $('<td><div id="arrow-'+arrowid+'">&#8592;</div>'+num.toString()+'</td>');
					$('#dpmatrix').append(cell);
					var idstring = 'r'+i +'-c'+ j ;
					cell.attr('id',idstring);

					if(i-1==j){
						//cell.attr('class',"diag");
					}
				}

			}
		}
		$('#dpmatrix').append("</tr>");

		//var cellstring = "#r"+(b.length+1)+"-c"+(a.length);
		//var lastcell = $(cellstring);
		//console.log(cellstring);
		//lastcell.attr('class',"optimal");
		checkpossible(alignment.Dir,b.length+1,a.length);

		// Trace_back
		for(i=a.length;a>2;i--){
			for(j=b.length;b>2;j--){

			}
		}
		if(a.length > 0 && b.length > 0){
			$('#result').empty();
			var seq_a_aligned = alignment.SeqA_aligned;
			var seq_b_aligned = alignment.SeqB_aligned;

			$('#result').append("Resulting Alignment : <br>");
			$('#result').append("<br>");
			$('#result').append(seq_a_aligned + " <br>");
			$('#result').append(getPattern(seq_a_aligned,seq_b_aligned) +" <br>");
			$('#result').append(seq_b_aligned + " <br>");
			$('#result').append("<br>");
			$('#result').append("Score : " + alignment.Score );

		}
		buildtable2();
	}
}

var delay = ( function() {
	var timer = 0;
	return function(callback, ms) {
		clearTimeout (timer);
		timer = setTimeout(callback, ms);
	};
})();

function checkpossible(tb,i,j){

	console.log(i + " , " + j);
	//if(i<=2 && j<=2) return;
	var cellstring = "#r"+i+"-c"+j;
	var lastcell = $(cellstring);
	//	lastcell.attr('class',"optimal",1000);
	if(animate){
		lastcell.addClass("optimal",(tb.length - i)*300);
	}else{
		lastcell.addClass("optimal");
	}
	console.log("value is : " + tb[j][i-1]);
	if(tb[j][i-1]===0){
		console.log("diag");
		checkpossible(tb,i-1,j-1);
	}else if(tb[j][i-1]==-1){
		console.log("left");
		checkpossible(tb,i,j-1);
	}else if (tb[j][i-1]==1) {
		console.log("up");
		checkpossible(tb,i-1,j);
	}
	else return;
}

function buildtable2(){
	// GCATGCU
	// GATTACA
	var a = $('input[name=seq1]').val();
	var b = $('input[name=seq2]').val();
	if (a.length > 0 && b.length >0){
		var alignment = Linear_alignment(a,b,1,-1,1);
		console.log(alignment);
		for (i = 0; i <= b.length+1; i++) {
			$('#dpmatrix2').append("<tr>");
			if(i<2){
				$('#dpmatrix2').append("<th>-</th>");
			}else{
				$('#dpmatrix2').append("<th>"+b.charAt(i-2)+"</th>");
			}
			for (j=0; j <= a.length; j++) {
				if(i===0){
					if(j===0){
						$('#dpmatrix2').append("<th>-</th>");
					}else{
						$('#dpmatrix2').append("<th>"+a.charAt(j-1)+"</th>");
					}
				}else{
					var num = alignment.Dir[j][i-1];

					var cell = $('<td>'+num.toString()+'</td>');
					$('#dpmatrix2').append(cell);
					var idstring = 'r'+i +"-c" + j ;
					cell.attr('id',idstring);
				}

			}
			$('#dpmatrix2').append("</tr>");

		}
	}
}


function getPattern(a_aligned,b_aligned){
	var str = "";
	for (var i = 0; i < a_aligned.length; i++) {
		if(a_aligned.charAt(i) == b_aligned.charAt(i)){
			str = str + "|";
		}else{
			str = str + ".";
		}
	}
	return str;
}

function randInt(){
	return Math.floor(Math.random() * 100) + 1;
}

function nextStep(){
	if(stepcounter < maxsteps){
		stepcounter  = stepcounter + 1 ;
		$("#counter").html("Step :" + stepcounter + " / " + maxsteps);
		highlight(stepcounter);
	}
	console.log(stepcounter);
}


function highlight(stepcounter){
	var modulo = stepcounter % 3;
	if(modulo === 0){
		$("#gap_a").addClass("highlight");
		$("#gap_b").removeClass("highlight");
		$("#match_missmatch").removeClass("highlight");
	}else if (modulo == 1) {
		$("#gap_b").addClass("highlight");
		$("#gap_a").removeClass("highlight");
		$("#match_missmatch").removeClass("highlight");
	} else if (modulo == 2) {
		$("#match_missmatch").addClass("highlight");
		$("#gap_a").removeClass("highlight");
		$("#gap_b").removeClass("highlight");
	}
}

function prevStep(){
	if(stepcounter > 0){
		stepcounter  = stepcounter - 1 ;
		$("#counter").html("Step :" + stepcounter + " / " + maxsteps);
		highlight(stepcounter);
	}
	console.log(stepcounter);
}

$(function () {
	var $slider = $('#slider');
	if ($slider.length > 0) {
		$slider.slider({
			min: 1,
			max: 5,
			value: 1,
			orientation: 'horizontal',
			range: 'min'
		}).addSliderSegments($slider.slider('option').max);
	}

	$("#slider").slider({
	  change: function( event, ui ) {
			gapcost = ui.value ;
			var gcost = $("#gcost");
			gcost.html(gapcost.toString());
			buildtable();
		}
	});

	$("select").select2({dropdownCssClass: 'dropdown-inverse'});


});

$(function () {
	var $slider2 = $('#slider2');
	if ($slider2.length > 0) {
		$slider2.slider({
			min: 1,
			max: 5,
			value: 1,
			orientation: 'horizontal',
			range: 'min'
		}).addSliderSegments($slider2.slider('option').max);
	}

	$("#slider2").slider({
		change: function( event, ui ) {
			scorematch = ui.value ;
			var smatch = $("#smatch");
			smatch.html(scorematch.toString());
			buildtable();
		}
	});

	$("select").select2({dropdownCssClass: 'dropdown-inverse'});


});

$(function () {
	var $slider3 = $('#slider3');
	if ($slider3.length > 0) {
		$slider3.slider({
			min: 0,
			max: 5,
			value: 1,
			orientation: 'horizontal',
			range: 'min'
		}).addSliderSegments($slider3.slider('option').max);
	}

	$("#slider3").slider({
		change: function( event, ui ) {
			scoremismatch = ui.value*(-1) ;
			var smismatch = $("#smismatch");
			smismatch.html(scoremismatch.toString());
			buildtable();
		}
	});

	$("select").select2({dropdownCssClass: 'dropdown-inverse'});


});

// Add segments to a slider
$.fn.addSliderSegments = function (amount, orientation) {
	return this.each(function () {
		if (orientation == "vertical") {
			var output = '', i;
			for (i = 1; i <= amount - 2; i++) {
				output += '<div class="ui-slider-segment" style="top:' + 100 / (amount - 1) * i + '%;"></div>';
			}
			$(this).prepend(output);
		} else {
			var segmentGap = 100 / (amount - 1) + "%", segment = '<div class="ui-slider-segment" style="margin-left: ' + segmentGap + ';"></div>';
			$(this).prepend(segment.repeat(amount - 2));
		}
	});
};
