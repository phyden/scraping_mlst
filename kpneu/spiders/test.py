# -*- coding: utf-8 -*-
import scrapy
from kpneu.items import KpneuItem

class TestSpider(scrapy.Spider):
    name = "test"
    allowed_domains = ["http://bigsdb.web.pasteur.fr"]
    start_urls = (
        'http://bigsdb.web.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef_public&page=downloadAlleles',
    )
    def parse(self, response):
        for sel in response.xpath('//div/div//td'):
            item = KpneuItem()
            item['title'] = sel.xpath('a/text()').extract()
            item['link'] = sel.xpath('a/@href').extract()
            item['desc'] = sel.xpath('text()').extract()
            if not len(item['link']) == 0 and len(item['desc']) == 0:
                url=item['link']
                url=str(url)[3:-2]
                url=response.urljoin(url)
                item["file_urls"] = [url]
                yield item
                print url



    def parse_dir_contents(self, response):
        print "test"
