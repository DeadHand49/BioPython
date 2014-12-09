

from __future__ import unicode_literals
import pubmedgrabber
import unittest


class PubMedTests(unittest.TestCase):

    def test_clean_titles_same(self):
        self.assertEqual(
            'Tauroursodeoxycholic Acid, a Bile Acid, Promotes Blood Vessel Repair by Recruiting Vasculogenic Progenitor Cells',
            pubmedgrabber.clean_title('Tauroursodeoxycholic Acid, a Bile Acid, Promotes Blood Vessel Repair by Recruiting Vasculogenic Progenitor Cells'))

    def test_clean_titles_in_vivo(self):
        self.assertEqual('A Stromal Cell-Free Culture System Generates Mouse Pro-T Cells that Can Reconstitute T-Cell Compartments <em>In Vivo</em>',
                         pubmedgrabber.clean_title('A stromal cell-free culture system generates mouse pro-T cells that can reconstitute T-cell compartments in vivo'))

    def test_clean_titles_phase_i(self):
        self.assertEqual('A Phase I Study of the Safety, Pharmacokinetics and Anti-Leukemic Activity of the Anti-CD123 Monoclonal Antibody CSL360 in Relapsed, Refractory or High-Risk Acute Myeloid Leukemia',
                         pubmedgrabber.clean_title('A Phase 1 study of the safety, pharmacokinetics and anti-leukemic activity of the anti-CD123 monoclonal antibody CSL360 in relapsed, refractory or high-risk acute myeloid leukemia'))

    def test_clean_titles_slash_cap(self):
        self.assertEqual('Myelodysplastic/Myeloproliferative',
                         pubmedgrabber.clean_title('Myelodysplastic/myeloproliferative'))

    def test_clean_titles_remove_britishness(self):
        self.assertEqual('Leukemia',
                         pubmedgrabber.clean_title('leukaemia'))

    def test_clean_titles_capital_after_colon(self):
        self.assertEqual('Engraftment: A',
                         pubmedgrabber.clean_title('Engraftment: a'))

    def test_clean_titles_in_vitro_dash(self):
        self.assertEqual('<em>In Vitro</em>-Induced',
                         pubmedgrabber.clean_title('in Vitro-Induced'))

    def test_clean_abstract_super(self):
        self.assertEqual('CD34<sup>+</sup>',
                         pubmedgrabber.clean_abstract('CD34+'))

    def test_clean_abstract_super_space(self):
        self.assertEqual('CD34<sup>+</sup>/Sca1<sup>+</sup> ',
                         pubmedgrabber.clean_abstract('CD34+ /Sca1+ '))

    def test_clean_abstract_comma_space(self):
        self.assertEqual('CD34<sup>+</sup>,',
                         pubmedgrabber.clean_abstract('CD34+ ,'))

    def test_clean_abstract_phase(self):
        self.assertEqual('Phase I',
                         pubmedgrabber.clean_abstract('phase 1'))

    def test_clean_abstract_number(self):
        self.assertEqual('nine',
                         pubmedgrabber.clean_abstract('9'))

    def test_clean_abstract_trailing(self):
        self.assertEqual('dose levels.',
                         pubmedgrabber.clean_abstract('dose levels. '))

    def test_clean_abstract_add_period(self):
        self.assertEqual('dose levels.',
                         pubmedgrabber.clean_abstract('dose levels '))

    def test_clean_abstract_remove_today(self):
        self.assertEqual(' ',
                         pubmedgrabber.clean_abstract(' today '))




if __name__ == '__main__':
    unittest.main()