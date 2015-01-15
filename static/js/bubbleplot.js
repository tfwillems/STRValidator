function draw_bubbleplot(inc_diag){
var h, halfh, halfw, margin, totalh, totalw, w;

h = 300;
w = 300;

margin = {
  left: 60,
  top: 40,
  right: 60,
  bottom: 40,
  inner: 5
};

halfh  = h + margin.top + margin.bottom;
totalh = halfh * 2;
halfw  = w + margin.left + margin.right;
totalw = halfw * 2;

function compute_scaling_factor(x, y, counts){
    var max_fac = 1000;
    var xdiff, ydiff, dist, fac;
    for (i = 0; i < x.length; i++){
	for (j = i+1; j < x.length; j++){
	    xdiff = x[i]-x[j];
	    ydiff = y[i]-y[j];
	    dist  = Math.sqrt(xdiff*xdiff + ydiff*ydiff);
	    fac   = dist/(Math.sqrt(counts[i]) + Math.sqrt(counts[j]));
	    max_fac = Math.min(fac, max_fac);
	}
    }
    return max_fac;
}

function gen_alignments(x, y) {
    $('#click_desc').html('');
    $('#alignments').html('<img src="/static/loading_icon.gif" style="display: block; margin: auto;">');
    var url  = "getalignments/" + x + "/" + y;
    $.get(url)
	.done(function(data) {
	    $("#alignments").html(data);
	})
	.fail(function(data) {
	      alert("Something went wrong" + data);
	});
}
 
$.get('/getbubbleinfo')
.done(function(data){
  var data = data.result.map(function(d){ return d3.values(d);});
  var min_x, max_x, min_y, max_y, min_v, max_v, chroms, starts, ends, samples;

  min_x   = d3.min(data.map(function(d){ return +d[0];}));
  max_x   = d3.max(data.map(function(d){ return +d[0];}));
  min_y   = d3.min(data.map(function(d){ return +d[1];}));
  max_y   = d3.max(data.map(function(d){ return +d[1];}));
  min_v   = Math.min(min_x, min_y);
  max_v   = Math.max(max_x, max_y);
  sizes   = data.map(function(d){ return +d[2];});
  bubblechart = scatterplot().xvar(0).yvar(1).xlab("lobSTR").ylab("Capillary").height(h).width(w).margin(margin)
      .pointcolor("orange").xNA({handle:false}).yNA({handle:false}).xlim([min_v-5,max_v+5]).ylim([min_v-5, max_v+5]).title("Bubble plot");
  d3.select("div#chart1").datum({
    data: data
  }).call(bubblechart);

  svg = d3.select("div#chart1").selectAll("svg");
  indtip = d3.tip().attr('class', 'd3-tip').html(function(d, i) {
     return "(" + d3.select(this).attr("x") + ", " + d3.select(this).attr("y") + "), ncalls = " + d3.select(this).attr("z");
  }).direction('e').offset([0, 10]);
  svg.call(indtip);
  
  return bubblechart.pointsSelect().on("mouseover", function(d) {
	  return d3.select(this).attr("fill", "red");
  }).on("mouseout", function(d) {
	  return d3.select(this).attr("fill", bubblechart.pointcolor());
  }).on("click", function(d) {
    gen_alignments(d3.select(this).attr("x"), d3.select(this).attr("y"));
    return 1;
  }).attr("x", function(d,i) {
    return data[i][0];
  }).attr("y", function(d,i) {
    return data[i][1];
  }).attr("z", function(d,i) {
    return sizes[i];
  }).attr("r", function(d,i){
    if (!inc_diag && d3.select(this).attr("x") == d3.select(this).attr("y"))
	 return 0;
     else
	 return scaling_fac*Math.sqrt(sizes[i]);
  }).on("mouseover.paneltip", indtip.show).on("mouseout.paneltip", indtip.hide);
})
.fail(function(data) {
   alert("Something went wrong" + data);
});

}
