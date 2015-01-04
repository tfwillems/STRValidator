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

import argparse
import os.path

# Set by command line options
bams      = ""
bais      = ""
fasta_dir = ""
vizalign  = ""
data_file = ""

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
    data = open("data/test_chrY.csv", "r")
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
    parser = argparse.ArgumentParser(description='Web server to compare STR calls and their associated alignments')
    parser.add_argument("--bams",     type=str, required=True, help='Comma-separated list of bam files')
    parser.add_argument("--bais",     type=str, required=True, help='Comma-separated list of bam index files. Order must match that of the --bams arguments')
    parser.add_argument("--fasta",    type=str, required=True, help='Directory containing chromosome fasta files')
    parser.add_argument("--vizalign", type=str, required=True, help='Full path for vizalign executable')
    parser.add_argument("--datafile", type=str, required=True, help="Full path for file containing bubble plot info in .csv format")
    args = parser.parse_args()
    
    # Check that .bam files exist
    bams = args.bams.split(",")
    for file in bams:
        if not os.path.isfile(file):
            exit("ERROR: BAM file %s does not exist"%(file))
    bams = ",".join(bams)

    # Check that .bai files exist
    bais = args.bais.split(",")
    for file in bais:
        if not os.path.isfile(file):
            exit("ERROR: BAI file %s does not exist"%(file))
    bais = ",".join(bais)

    # Check that fasta directory exists
    if not os.path.isdir(args.fasta):
        exit("ERROR: Fasta directory does not exist")
    fasta_dir = args.fasta

    # Check that vizalign path exists
    if not os.path.isfile(args.vizalign):
        exit("ERROR: Invalid vizalign path")
    vizalign = args.vizalign

    # Check that the data file exists
    if not os.path.isfile(args.datafile):
        exit("ERROR: Invalid data file path")
    data_file = args.datafile

    # Run flask app
    app.run(host='0.0.0.0', port=6015, debug=False)
