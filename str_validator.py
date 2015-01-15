#!/usr/bin/env python2.7

from flask import Flask, url_for, make_response, render_template, request
from flask import json, jsonify, current_app, Blueprint
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy.sql import func

from json import dumps
from sqlalchemy.orm import class_mapper
from sqlalchemy import or_, and_, ForeignKey
from sqlalchemy import func

from flask.ext.compress import Compress
from flask.ext.assets import Environment, Bundle

from subprocess import Popen,PIPE

import argparse
import os.path
import collections

admin = Blueprint('admin', __name__, static_folder='static', url_prefix='/')
db    = SQLAlchemy()

class HaploidCalls(db.Model):
    __tablename__ = 'haploid_calls'
    id         = db.Column(db.Integer, primary_key=True)
    chrom      = db.Column(db.Text)
    start      = db.Column(db.Integer)
    end        = db.Column(db.Integer)
    total_one  = db.Column(db.Integer)
    total_two  = db.Column(db.Integer)
    g1_1       = db.Column(db.Integer)
    g1_2       = db.Column(db.Integer)
    g2_1       = db.Column(db.Integer)
    g2_2       = db.Column(db.Integer)
    sample_1   = db.Column(db.Text)
    sample_2   = db.Column(db.Text)

class DiploidCalls(db.Model):
    __tablename__ = 'diploid_calls'
    id         = db.Column(db.Integer, primary_key=True)
    chrom      = db.Column(db.Text)
    start      = db.Column(db.Integer)
    end        = db.Column(db.Integer)
    total_one  = db.Column(db.Integer)
    total_two  = db.Column(db.Integer)
    g1_1       = db.Column(db.Integer)
    g1_2       = db.Column(db.Integer)
    g2_1       = db.Column(db.Integer)
    g2_2       = db.Column(db.Integer)
    sample_1   = db.Column(db.Text)
    sample_2   = db.Column(db.Text)

@admin.route('/')
def index():
    return render_template("new_index.html")

def get_alignment_info(chroms, starts, stops, sample_lists, annot_lists):
    input        = "\n".join(map(lambda x: "\t".join(x), zip(chroms, starts, stops, sample_lists, annot_lists)))
    cmd          = "%s --bams %s --indexes %s --out /dev/stdout --regions /dev/stdin --fasta %s"%(current_app.config['vizalign'], current_app.config['bams'],
                                                                                                  current_app.config['bais'], current_app.config['fasta_dir'])
    stream       = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE)
    res          = stream.communicate(input=input)[0]
    return res

@admin.route('getbubbleinfo')
def get_bubble_info():
    # Compute bubble counts if not precomputed
    if current_app.config['BUBBLE_COUNTS'] is None:
        hap_counts = collections.defaultdict(int)
        for record in db.session().query(HaploidCalls):
            hap_counts[(record.total_one, record.total_two)] += 1
        current_app.config['BUBBLE_COUNTS'] = hap_counts        
    vals = map(lambda x: [x[0][0], x[0][1], x[1]], current_app.config['BUBBLE_COUNTS'].items())
    return jsonify(result=vals)

@admin.route('getalignments/<x>/<y>')
def get_alignments(x, y):
    x, y = int(x), int(y)
    res = db.session.query(HaploidCalls).filter(HaploidCalls.total_one == x, HaploidCalls.total_two == y)
    chroms, starts, stops, sample_lists, annot_lists = [], [], [], [], []
    for record in res:
        chroms.append(record.chrom)
        starts.append(str(record.start))
        stops.append(str(record.end))
        if record.sample_1 == record.sample_2:
            sample_lists.append(record.sample_1)
            annot_lists.append(str(record.g1_1)+"/"+str(record.g1_2))
        else:
            sample_lists.append(record.sample_1 + "," + record.sample_2)
            annot_lists.append(str(record.g1_1)+"/"+str(record.g1_2)+ "," + str(record.g2_1)+"/"+str(record.g2_2))

    res      = get_alignment_info(chroms, starts, stops, sample_lists, annot_lists)
    response = make_response(res)
    response.headers.add('Access-Control-Allow-Origin','*')
    return response


def create_app(bams, bais, fasta_dir, vizalign, data_file):
    app = Flask(__name__)

    # Check that .bam files exist
    bams = bams.split(",")
    for file in bams:
        if not os.path.isfile(file):
            exit("ERROR: BAM file %s does not exist"%(file))
    app.config['bams'] = ",".join(bams)

    # Check that .bai files exist
    bais = bais.split(",")
    for file in bais:
        if not os.path.isfile(file):
            exit("ERROR: BAI file %s does not exist"%(file))
    app.config['bais'] = ",".join(bais)

    # Check that fasta directory exists
    if not os.path.isdir(fasta_dir):
        exit("ERROR: Fasta directory does not exist")
    app.config['fasta_dir'] = fasta_dir

    # Check that vizalign path exists
    if not os.path.isfile(vizalign):
        exit("ERROR: Invalid vizalign path")
    app.config['vizalign'] = vizalign

    # Check that the data file exists
    if not os.path.isfile(data_file):
        exit("ERROR: Invalid file path for bubble plot data")
    app.config['data_file'] = data_file

    app.register_blueprint(admin)
    Compress(app) # Compresses large responses (useful for alignment viewing)
    assets = Environment(app)
    js     = Bundle('js/d3.min.js',    
                    'js/colorbrewer.js',
                    'js/d3-tip.min.js',  
                    'js/panelutil.js',  
                    'js/scatterplot.js',
                    'js/jquery-2.1.3.min.js',
                    'js/bubbleplot.js',
                    filters='jsmin', output='gen/packed.js')
    assets.register('str_validator_js_all', js)

    app.config['BUBBLE_COUNTS']           = None
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///call_data.sqlite'
    db.init_app(app)
    with app.app_context():
        db.create_all()

    return app


if __name__== '__main__':
    parser = argparse.ArgumentParser(description='Web server to compare STR calls and their associated alignments')
    parser.add_argument("--bams",     type=str, required=True, help='Comma-separated list of bam files')
    parser.add_argument("--bais",     type=str, required=True, help='Comma-separated list of bam index files. Order must match that of the --bams arguments')
    parser.add_argument("--fasta",    type=str, required=True, help='Directory containing chromosome fasta files')
    parser.add_argument("--vizalign", type=str, required=True, help='Full path for vizalign executable')
    parser.add_argument("--datafile", type=str, required=True, help="Full path for file containing bubble plot info in .csv format")
    args = parser.parse_args()

    # Create app
    my_app = create_app(args.bams, args.bais, args.fasta, args.vizalign, args.datafile)

    # Run it
    my_app.run(host='0.0.0.0', port=6015, debug=True)
