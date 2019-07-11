#!/usr/bin/env python

import scrapy
import os
import re
from io import StringIO
from Bio import SeqIO

outputdir="lmono_mlst"

def check_sequences(response):
        """
        Function to check for sanity:
            no two ids in the same file may be identical
            no two sequences in the same file may be identical
            sequences must be nucleotide only!
        """
        ids_seen = set()
        sequences_seen = set()
        write_sequences = []
        fasta_io = StringIO(response.body.decode('ascii'))
        for record in SeqIO.parse(fasta_io,"fasta"):
                write = True
                if record.id in ids_seen:
                        write = False
                if str(record.seq) in sequences_seen:
                        write = False
                ids_seen.add(record.id)
                sequences_seen.add(str(record.seq))
                if write:
                        write_sequences.append(record)

        for record in write_sequences:
                if re.search("[^ACGT]",str(record.seq)):
                        write_sequences.drop(record)
                        print("INFO: %s contains forbidden character, dropping sequence" % record.id)

        fasta_io.close()

        return write_sequences

class SchemeSpider(scrapy.Spider):
        name = "scheme"
        allowed_domains = ['bigsdb.pasteur.fr']
        def start_requests(self):
                urls = [
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=5&render=1', # Virulence
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=6&render=1', # AB resistance
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=7&render=1', # Metal & Detergent resistance
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=8&render=1', # Stress Islands
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=10&render=1', # Listeria Stress Islands
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=11&render=1', # sigB operon
                        'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef&page=downloadAlleles&scheme_id=12&render=1' # Rhamnose operon
                ]
                for url in urls:
                        yield scrapy.Request(url=url, callback=self.parse)

        def parse(self, response):
                page = response.url.split("/")[:3]
                scheme_name = response.xpath('//h2/a/text()').get()
                dirname = scheme_name.replace(" ","_").replace("&","")
                if not os.path.exists(outputdir+"/"+dirname):
                        os.makedirs(outputdir+"/"+dirname)
                for gene in response.xpath("//table[@class='resultstable']/tr/td/a[@href]"):
                        link = gene.attrib["href"]
                        if "download" in link:
                                fasta_link = "/".join(page) + link
                                print(fasta_link)
                                yield scrapy.Request(url=fasta_link, callback=self.store_fasta, meta={"dirname": dirname})


        def store_fasta(self, response):
                dirname=response.meta["dirname"]
                gene_name = response.url.split("=")[-1]
                write_sequences = check_sequences(response)
                SeqIO.write(write_sequences, outputdir+"/"+dirname+"/"+gene_name+".fasta","fasta")
