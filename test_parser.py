import unittest
from parser import Parser


class TestParser(unittest.TestCase):
    def setUp(self):
        self.parser = Parser()
        self.parser.read_ecs()
        self.parser.read_embeddings()
        self.assertTrue(len(self.parser.uni_proteins) > 0)

    def tearDown(self):
        self.parser = None

    def test_find_by_entry(self):
        protein = self.parser.find_by_entry('A0A009IHW8')
        self.assertIsNotNone(protein)

    def test_find_all(self):
        proteins = self.parser.find_all()
        self.assertIsNotNone(proteins)


if __name__ == '__main__':
    unittest.main()
