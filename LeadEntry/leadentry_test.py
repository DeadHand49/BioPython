# coding=utf-8
from __future__ import unicode_literals

__author__ = 'memery'

import leadentry
import unittest
import urllib2


class LeadEntryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.url = 'http://www.mesenchymalcellnews.com/issue/volume-6-43-nov-18/'
        cls.newsletter = leadentry.Newsletter(cls.url)
        cls.article = cls.newsletter.articles[0]
        cls.author = cls.article.authors[-1]

    def test_parse_connexon_proper_length(self):
        self.assertEqual(12, len(self.newsletter.articles))

    def test_newsletter_specific_lead_source(self):
        self.assertEqual('Mesenchymal Cell News 6.43', self.newsletter.info['Specific Lead Source'])

    def test_newsletter_archived_link(self):
        self.assertEqual('mesenchymalcellnews.com/issue/volume-6-43-nov-18/',
                         self.newsletter.info['Newsletter Archived Link'])

    def test_newsletter_lead_source(self):
        self.assertEqual('Connexon', self.newsletter.info['Lead Source'])

    def test_get_publication_date(self):
        self.assertEqual('11/10/2014', self.article.info['Publication Date'])

    def test_get_publication_date_no_day(self):
        self.assertEqual('11/29/2014', self.newsletter.articles[8].info['Publication Date'])

    def test_publication_link(self):
        test_open = urllib2.urlopen('http://www.pnas.org/content/early/2014/11/06/1416121111.abstract')
        self.assertIn(self.article.info['Publication Link'], test_open.url)

    def test_publication_link_no_resolve(self):
        self.assertEqual('DOI cannot be resolved: http://dx.doi.org/10.1002/stem.1902',
                         self.newsletter.articles[2].info['Publication Link'])

    def test_article_title(self):
        self.assertEqual('TSG-6 as a Biomarker to Predict Efficacy of Human Mesenchymal '
                         'Stem/Progenitor Cells (hMSCs) in Modulating Sterile Inflammation In Vivo'.lower(),
                         self.article.info['Article Title'].lower())

    def test_author_first_name(self):
        self.assertEqual('Darwin', self.author.info['First Name'])

    def test_author_last_name(self):
        self.assertEqual('Prockop', self.author.info['Last Name'])

    def test_author_email(self):
        self.assertEqual('Prockop@medicine.tamhsc.edu', self.author.info['Email'])

    def test_author_company(self):
        self.assertIn('Texas A&M Health Science Center', self.author.info['Company'])

    def test_author_department(self):
        self.assertEqual('Institute for Regenerative Medicine', self.author.info['Department'])