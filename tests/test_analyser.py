import os
import unittest
from analyser import Analyser
from parser import Parser


class TestParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Shenanigans to get the path to the data files to work when running via unit tests
        cls.script_directory = os.path.dirname(os.path.abspath(__file__))
        class_under_test_directory = os.path.join(cls.script_directory, '..')
        os.chdir(class_under_test_directory)

        cls.parser = Parser();
        cls.parser.parse()
        cls.analyser = Analyser()

    @classmethod
    def tearDownClass(cls):
        cls.analyser = None
        cls.parser = None

    def test_analyse(self):
        #  def analyse(self, uniproteins, file_prefix):
        #  test first 10
        if len(self.parser.uni_proteins) > 0:
            first_10 = list(self.parser.uni_proteins.values())[:10]
            self.analyser.analyse(first_10, 'first_ten')
        else:
            print("No uni_proteins were parsed")


if __name__ == '__main__':
    unittest.main()
