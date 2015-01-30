# coding=utf-8
"""Unittests for leadentry.py"""

from __future__ import unicode_literals

__author__ = 'memery'

import leadentry
import unittest
import mock
from bs4 import BeautifulSoup, element


class LeadEntryTest(unittest.TestCase):

    def setUp(self):
        self.batch = leadentry.Batch('matthew.emery@stemcell.com')
        self.batch.pmids = ['24411336']
        self.batch.pubmed_xml = self.batch.create_pubmed_xml()
        self.batch.parse_pubmed_soup()

    def test_batch_title(self):
        self.assertEqual('A defined xeno-free and feeder-free culture system for the derivation, '
                         'expansion and direct differentiation of transgene-free patient-specific '
                         'induced pluripotent stem cells', self.batch.articles[0].info['Article Title'])

    def test_batch_abstract(self):
        self.assertIn('A defined xeno-free system for patient-specific iPSC derivation and differentiation',
                      self.batch.articles[0].info['Abstract'])

    def test_batch_pub_link(self):
        self.assertEqual('http://www.sciencedirect.com/science/article/pii/S0142961213015342',
                         self.batch.articles[0].info['Publication Link'])


class NewsletterTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tester = leadentry.Newsletter('http://www.mesenchymalcellnews.com/issue/volume-7-03-jan-27/',
                                          'matthew.emery@stemcell.com')

    def test_get_specific_lead_source(self):
        self.assertEqual('Mesenchymal Cell News 7.03', self.tester.find_specific_lead_source())

    def test_get_search_term(self):
        self.assertEqual('Connexon; Mesenchymal Cell News', self.tester.info['Search Term'])


class ArticleTest(unittest.TestCase):
    def setUp(self):
        tester = leadentry.Batch('matthew.emery@stemcell.com')
        with open('leadentry_test.xml') as xml:
            tester.pubmed_xml = BeautifulSoup(xml.read())
            self.test_article = leadentry.Article(info={'PMID': '25433608',
                                                        'Article Title': 'Mesenchymal stromal cells form vascular '
                                                                         'tubes when placed in fibrin sealant and '
                                                                         'accelerate wound healing in vivo'})
            tester.add_article(self.test_article)
            tester.parse_pubmed_soup()

    def test_in_info(self):
        self.assertTrue(self.test_article.in_info('PMID'))

    def test_not_in_info(self):
        self.assertFalse(self.test_article.in_info('Cell Type'))

    def test_get_info(self):
        self.assertEqual(self.test_article.get_info('PMID'),
                         '25433608')

    def test_no_get_info(self):
        self.assertIsNone(self.test_article.get_info('Cell Type'))

    def test_look_up_pmid(self):
        self.assertEqual('25433608',
                         self.test_article.lookup_up_pmid(self.test_article.get_info('Article Title'),
                                                          'matthew.emery@stemcell.com'))

    def test_look_up_pmid_first_fail(self):
        """Need to find an example of a test that fails at first, but succeeds on Britishness"""
        pass

    def test_look_up_pmid_complete_fail(self):
        """Need to find an example that fails in completely and returns None"""
        pass

    def test_translate_british(self):
        pub_example = 'Leukemia Tumor Signaling α'
        self.assertEqual('Leukaemia Tumour Signalling alpha',
                         self.test_article.translate_british(pub_example))

    def test_update_info_dict(self):
        self.test_article.update_info_dict('PMID', 'test')
        self.assertDictContainsSubset('PMID', 'test')

    def test_get_info_items(self):
        self.assertEqual((('PMID', '25433608'),
                         ('Article Title', 'Mesenchymal stromal cells form vascular tubes when placed in fibrin '
                                           'sealant and accelerate wound healing in vivo')),
                         self.test_article.get_info_items())

    def test_find_date(self):
        self.assertEqual('11/26/2014', self.test_article.find_date())

    def test_cant_find_date(self):
        """Need to find example where find date fails"""
        pass

    def test_find_pmid_in_xml(self):
        self.assertIsInstance(self.test_article.info['Tag'], element.Tag)


class AuthorTest(unittest.TestCase):
    pass

if __name__ == "__main__":
    unittest.main()