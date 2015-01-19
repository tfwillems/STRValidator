
function new_draw_barplot(inc_match, min_val, max_val) {
    function gen_alignments(diff) {
	$('#click_desc').html('');
	$('#alignments').html('<img src="/static/loading_icon.gif" style="display: block; margin: auto;">');
	var url  = "getdiffalignments/" + diff;
	$.get(url)
	    .done(function(data) {
		    $("#alignments").html(data);
		})
	    .fail(function(data) {
		    alert("Something went wrong" + data);
		});
    }

    var chart;
    nv.addGraph(function() {
	    chart = nv.models.multiBarChart()
		.transitionDuration(150)
		.reduceXTicks(false)   //If 'false', every single x-axis tick label will be rendered.
		.rotateLabels(0)      //Angle to rotate x-axis labels.
		.showControls(false)  //Allow user to switch between 'Grouped' and 'Stacked' mode.
		.groupSpacing(0.1)    //Distance between each group of bars.
		.showLegend(false)
		.tooltips(true)
		.tooltip( function(key, x, y, e, graph) { return '<p>' + "# calls = " + y + '</p>'});

	    chart.xAxis
		.tickFormat(d3.format(',d'))
		.axisLabel("Difference (bp)");
	    chart.yAxis
		.tickFormat(d3.format(',d'))
		.axisLabel("Number of calls");

	    $.get('/getdiffinfo')
		.done(function(data){
			var data    = data.result;
			var gt      = d3.sum(data.filter(function(x) { return +x[0] >= max_val; }), function(x) { return +x[1];});
			var lt      = d3.sum(data.filter(function(x) { return +x[0] <= min_val; }), function(x) { return +x[1];});
			data = data.filter(function(x) { return (inc_match || +x[0] != 0) && +x[0] > min_val && +x[0] < max_val;});
			if(lt != 0)
			    data.push([min_val, lt]);
		        if(gt != 0)
			    data.push([max_val, gt]);
			data     = data.sort(function(a,b){ return +a[0] - +b[0]; });
			data     = data.map(function(d){ return {"x": +d[0], "y": +d[1]};  });
			data     = [{key:"Dataset", values: data}];
			d3.select('#bar_plot_container svg')
			    .datum(data)
			    .call(chart);
			nv.utils.windowResize(chart.update);
		    });
	    
	    chart.multibar.dispatch.on('elementClick', function(e){ gen_alignments(e.point.x); });
	    return chart;
	});
}