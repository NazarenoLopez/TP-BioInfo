import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from Ex1 import procesar_secuencia_y_detectar_orf  # noqa: E402

GBK_INPUT = ROOT / "HTT_mRNA.gbk"


class Ex1Test(unittest.TestCase):
    """Tests del Ej1: generación de los 6 marcos y detección del frame correcto."""

    def setUp(self) -> None:
        self.assertTrue(
            GBK_INPUT.exists(),
            f"No se encuentra el archivo de entrada {GBK_INPUT}",
        )

    def test_genera_seis_marcos(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.fas"
            with contextlib.redirect_stdout(io.StringIO()):
                procesar_secuencia_y_detectar_orf(str(GBK_INPUT), str(out))
            with open(out) as handle:
                records = list(SeqIO.parse(handle, "fasta"))

        self.assertEqual(len(records), 6, "Deberían generarse exactamente 6 marcos")
        ids_esperados = {
            f"Frame_{i}_{sentido}"
            for i in (1, 2, 3)
            for sentido in ("Fwd", "Rev")
        }
        self.assertEqual({r.id for r in records}, ids_esperados)

    def test_detecta_frame_correcto_HTT(self) -> None:
        # En HTT_mRNA.gbk el CDS arranca en posición 146 (146 mod 3 = 2),
        # por lo tanto el marco correcto debe ser Frame_2_Fwd.
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.fas"
            with contextlib.redirect_stdout(io.StringIO()):
                marco = procesar_secuencia_y_detectar_orf(str(GBK_INPUT), str(out))
        self.assertEqual(marco, "Frame_2_Fwd")


if __name__ == "__main__":
    unittest.main()
