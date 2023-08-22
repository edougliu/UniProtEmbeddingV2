import unittest
from parser import Parser


class TestParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parser = Parser()
        cls.parser.read_ecs()
        cls.parser.read_embeddings()

    @classmethod
    def tearDownClass(cls):
        cls.parser = None

    def test_find_by_entry(self):
        protein1 = self.parser.find_by_entry('A0A009IHW8')
        self.assertIsNotNone(protein1)
        self.assertEqual(protein1.entry, 'A0A009IHW8')
        self.assertIsNotNone(protein1.embedding)

    def test_find_all(self):
        proteins = self.parser.find_all()
        self.assertIsNotNone(proteins)
        self.assertEqual(len(proteins), 569793)


if __name__ == '__main__':
    unittest.main()
