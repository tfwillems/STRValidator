#!/bin/bash

# Generate input files chrY_stats and autosomal_stats

# Construct SQL database
SQL="
CREATE TABLE haploid_calls (
id integer primary key,
chrom TEXT,
start NUMERIC,
end NUMERIC,
total_one NUMERIC,
total_two NUMERIC,
g1_1 NUMERIC,
g1_2 NUMERIC,
g2_1 NUMERIC,
g2_2 NUMERIC,
sample_1 TEXT,
sample_2 TEXT
);

CREATE TABLE diploid_calls (
id integer primary key,
chrom TEXT,
start NUMERIC,
end NUMERIC,
total_one NUMERIC,
total_two NUMERIC,
g1_1 NUMERIC,
g1_2 NUERMIC,
g2_1 NUMERIC,
g2_2 NUMERIC,
sample_1 TEXT,
sample_2 TEXT 
);

.separator \"	\"
.import chrY_stats haploid_calls
CREATE INDEX hap_index1 ON haploid_calls(total_one, total_two);
CREATE INDEX hap_index2 ON haploid_calls(chrom ASC);

.import autosomal_stats diploid_calls
CREATE INDEX dip_index1 ON diploid_calls(total_one, total_two);
CREATE INDEX dip_index2 ON diploid_calls(chrom ASC);
"
rm -f call_data.sqlite
echo "$SQL" | sqlite3 call_data.sqlite
