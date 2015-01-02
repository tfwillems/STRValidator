#!/usr/bin/env python2.7

from flask import Flask, url_for, make_response, render_template, request
from flask import Flask, json, jsonify
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy.sql import func

from json import dumps
from sqlalchemy.orm import class_mapper
from sqlalchemy import or_, and_, ForeignKey

from flask.ext.compress import Compress
from flask.ext.assets import Environment, Bundle

from subprocess import Popen,PIPE

app = Flask(__name__)

# Compresses large responses (useful for alignment viewing)
Compress(app)

assets = Environment(app)
js     = Bundle('js/d3.min.js',    
                'js/colorbrewer.js',
                'js/d3-tip.min.js',  
                'js/panelutil.js',  
                'js/scatterplot.js',
                'js/jquery-2.1.3.min.js',
                filters='jsmin', output='gen/packed.js')
assets.register('js_all', js)

bam_dir    = "/Users/tfwillems/Desktop/Club_files/chry/final_callset/chrY_BAMs/filtered_bams/"
fasta_dir  = "/Users/tfwillems/Desktop/Coding/dbase/"
bam_suffix = ".chrY.bam"
bai_suffix = ".chrY.bam.bai"
vizalign   = "/Users/tfwillems/Desktop/Coding/HipSTR/vizalign"
pops = ["ASW", "CEU", "CHB", "CHS", "CLM", "FIN", "GBR", "IBS", "JPT", "LWK", "MXL", "PUR", "TSI", "YRI"]
pops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", 
        "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
bams = ",".join(map(lambda x: bam_dir + x + bam_suffix, pops))
bais = ",".join(map(lambda x: bam_dir + x + bai_suffix, pops))


@app.route('/')
def index():
    return render_template("new_index.html")

def get_alignment_info(chroms, starts, stops, sample_lists):
    chroms       = chroms.split(",")
    starts       = starts.split(",")
    stops        = stops.split(",")
    sample_lists = map(lambda x: x.replace("-",","), sample_lists.split(","))
    input        = "\n".join(map(lambda x: "\t".join(x), zip(chroms, starts, stops, sample_lists)))
    print(input)
    cmd    ="%s --bams %s --indexes %s --out /dev/stdout --regions /dev/stdin --fasta %s"%(vizalign, bams, bais, fasta_dir)
    stream = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE)
    res    = stream.communicate(input=input)[0]
    return res

@app.route('/getbubbleinfo')
def get_bubble_info():
    info = ""
    data = open("test_chrY.csv", "r")
    for line in data:
        info += line
    data.close()
    return info

@app.route('/getalignments/<chroms>/<starts>/<stops>/<sample_lists>')
def get_alignments(chroms, starts, stops, sample_lists):
    res = get_alignment_info(chroms, starts, stops, sample_lists)
    response = make_response(res)
    response.headers.add('Access-Control-Allow-Origin','*')                                                                                                                      
    return response

if __name__== '__main__':
    app.run(host='0.0.0.0', port=6015, debug=False)
