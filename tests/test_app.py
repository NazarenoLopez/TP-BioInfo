import unittest

from src.app import get_message


class AppTest(unittest.TestCase):
    def test_get_message(self) -> None:
        self.assertEqual(get_message(), "Hola Mundo desde App")


if __name__ == "__main__":
    unittest.main()
