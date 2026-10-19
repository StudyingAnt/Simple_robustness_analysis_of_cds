"""Check complete coverage and sequential execution with fewer submitted jobs."""
import argparse
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '06_Bat_genome'))
import run_bat_cosmic_tasks as tasks


class BatTaskPlanTests(unittest.TestCase):
    def test_fifty_jobs_cover_all_two_hundred_original_cases_once(self):
        plan = [tasks.task_cases(i, 50, 100) for i in range(50)]
        self.assertTrue(all(len(cases) == 4 for cases in plan))
        self.assertEqual(plan[0], [0, 1, 100, 101])
        self.assertEqual(plan[-1], [98, 99, 198, 199])
        flattened = [i for cases in plan for i in cases]
        self.assertEqual(sorted(flattened), list(range(200)))
        for cases in plan:
            self.assertEqual([tasks.array_case(i, 100)[0] for i in cases],
                             ['molMol2', 'myoMyo6', 'molMol2', 'myoMyo6'])

    def test_uneven_task_count_keeps_full_coverage_and_stable_retry_mapping(self):
        plan = [tasks.task_cases(i, 7, 100) for i in range(7)]
        self.assertEqual(sorted(i for cases in plan for i in cases), list(range(200)))
        self.assertEqual(tasks.task_cases(3, 7, 100), plan[3])

    def test_invalid_indices_and_counts_fail(self):
        for values in [(-1, 50, 100), (50, 50, 100), (0, 0, 100), (0, 101, 100), (0, 1, 0)]:
            with self.assertRaises(ValueError):
                tasks.task_cases(*values)

    def test_failure_stops_later_cases_and_retry_uses_the_same_cases_and_settings(self):
        args = argparse.Namespace(task_index=0, task_count=2, shards_per_species=4,
                                  n_traj=1000, seed=42, workers=96)
        visited = []

        def fail_second(case):
            visited.append(case.array_index)
            if len(visited) == 2:
                raise RuntimeError('interrupted')

        with patch.object(tasks, 'run', side_effect=fail_second), redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError, 'interrupted'):
                tasks.run_task(args)
        self.assertEqual(visited, [0, 1])
        retried = []
        with patch.object(tasks, 'run', side_effect=lambda case: retried.append(
                (case.array_index, case.n_traj, case.seed, case.workers))), redirect_stdout(io.StringIO()):
            tasks.run_task(args)
        self.assertEqual(retried, [(i, 1000, 42, 96) for i in (0, 1, 4, 5)])
        self.assertFalse(hasattr(args, 'array_index'))


if __name__ == '__main__':
    unittest.main()
