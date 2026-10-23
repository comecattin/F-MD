"""Configuration contract: defaults, comments, and invalid physical inputs."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('input_parser', Path(__file__).resolve().parents[1] / 'src/input_parser.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def test_defaults(self):
        self.assertEqual(parser.config_values({}), (30, 1000, 0.0001, 10.0, 1))

    def test_comments_and_whitespace(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'input.in'
            source.write_text('# comment\n\n n_atoms : 6 # two waters\ndt: 1e-5\nstride: 2\n')
            self.assertEqual(parser.config_values(parser.parse_config(source)), (6, 1000, 1e-5, 10.0, 2))

    def test_invalid_values(self):
        for key, values in {'n_atoms': ['0', '-3', '4', '3.5'], 'n_steps': ['0', '-1'],
                            'stride': ['0', '-1'], 'dt': ['0', '-0.1', 'nan', 'inf'],
                            'box_length': ['0', '-1', 'nan', 'inf']}.items():
            for value in values:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    parser.config_values({key: value})

    def test_malformed_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'input.in'
            source.write_text('n_atoms 6\n')
            with self.assertRaisesRegex(ValueError, 'Line 1'):
                parser.parse_config(source)


if __name__ == '__main__':
    unittest.main()
