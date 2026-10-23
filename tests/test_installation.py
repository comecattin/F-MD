"""Exercise the public commands after installation, outside the source tree."""
import importlib.util
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


class InstallationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.installation = tempfile.TemporaryDirectory(prefix='fmd installation ')
        cls.addClassCleanup(cls.installation.cleanup)
        cls.prefix = Path(cls.installation.name) / "prefix with 'quotes'"
        subprocess.run(['make', '-j2', 'install', f'PREFIX={cls.prefix}'], cwd=ROOT, check=True)
        cls.env = dict(os.environ, FMD_PYTHON=sys.executable, MPLBACKEND='Agg')
        cls.env['PATH'] = str(cls.prefix / 'bin') + os.pathsep + cls.env['PATH']

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(prefix='fmd simulation ')
        self.addCleanup(self.work.cleanup)
        self.directory = Path(self.work.name)
        self.env = dict(self.env, MPLCONFIGDIR=str(self.directory / 'mpl'))

    def run_command(self, *args, check=True):
        return subprocess.run(args, cwd=self.directory, env=self.env, capture_output=True,
                              text=True, check=check, timeout=60)

    def simulate(self):
        config = self.directory / "input with 'quotes'.in"
        config.write_text('n_atoms: 3\nn_steps: 40\ndt: 0.00001\nbox_length: 20\nstride: 2\n')
        result = self.run_command('md_simulation', str(config))
        self.assertIn('Simulation finished', result.stdout)

    def test_installed_simulation_outputs(self):
        self.simulate()
        lines = (self.directory / 'energies.dat').read_text().splitlines()
        self.assertEqual(len(lines), 20 * 4)
        energies = []
        for index in range(20):
            frame = lines[index * 4:(index + 1) * 4]
            self.assertEqual(int(frame[0].split(':')[1]), 2 * (index + 1))
            values = [float(line.split(':')[1]) for line in frame[1:]]
            self.assertTrue(all(math.isfinite(value) for value in values))
            self.assertAlmostEqual(values[0] + values[1], values[2], places=10)
            energies.append(values[2])
        self.assertLess(max(energies) - min(energies), 1.0e-3)
        lines = (self.directory / 'trajectories.dat').read_text().splitlines()
        self.assertEqual(len(lines), 21 * 4)  # Includes the initial frame.
        for index, line in enumerate(lines):
            if index % 4 == 0:
                self.assertEqual(int(line.split(':')[1]), (index // 4) * 2)
            else:
                xyz = list(map(float, line.split()))
                self.assertEqual(len(xyz), 3)
                self.assertTrue(all(math.isfinite(value) for value in xyz))
        self.assertEqual(len((self.directory / 'trajectories.arc').read_text().splitlines()), 21 * 4)
        self.assertFalse((self.directory / 'parsed_config.tmp').exists())

    def test_missing_input_and_bad_arguments_fail(self):
        self.assertNotEqual(self.run_command('md_simulation', check=False).returncode, 0)
        self.assertNotEqual(self.run_command('md_simulation', 'missing.in', check=False).returncode, 0)
        # A stale parser output must not allow invalid input to launch a simulation.
        (self.directory / 'parsed_config.tmp').write_text('6\n40\n0.00001\n20\n2\n')
        config = self.directory / 'bad.in'
        config.write_text('stride: 0\n')
        self.assertNotEqual(self.run_command('md_simulation', str(config), check=False).returncode, 0)
        self.assertFalse((self.directory / 'energies.dat').exists())

    def test_installed_energy_plot(self):
        self.simulate()
        self.run_command('fmd-plot-energies', 'energies.dat')
        self.assertTrue((self.directory / 'energies.pdf').read_bytes().startswith(b'%PDF'))

    def test_trajectory_gif(self):
        self.simulate()
        # Test GIF rendering even on development machines without FFmpeg.
        spec = importlib.util.spec_from_file_location('trajectories', ROOT / 'src/plot_trajectories.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.steps, module.positions = module.read_trajectories(self.directory / 'trajectories.dat')
        module.x_max = module.y_max = 20.0
        fig, module.ax = module.plt.subplots()
        self.addCleanup(module.plt.close, fig)
        animation = module.FuncAnimation(fig, module.animate, frames=2, repeat=False)
        path = self.directory / 'trajectory.gif'
        animation.save(path, writer=module.PillowWriter(fps=24))
        with Image.open(path) as image:
            self.assertEqual(image.format, 'GIF')

    def test_installed_trajectory_exports(self):
        if not shutil.which('ffmpeg'):
            if os.environ.get('FMD_REQUIRE_FFMPEG') == '1':
                self.fail('FFmpeg is required for the CI installation test')
            self.skipTest('Install FFmpeg to test the installed GIF/MP4 command')
        self.simulate()
        self.run_command('fmd-plot-trajectories')
        with Image.open(self.directory / 'particle_trajectories.gif') as image:
            self.assertEqual(image.format, 'GIF')
        video = self.directory / 'particle_trajectories.mp4'
        self.assertGreater(video.stat().st_size, 0)
        self.assertIn(b'ftyp', video.read_bytes()[:32])


if __name__ == '__main__':
    unittest.main()
