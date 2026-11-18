"""Helper for test_dag.py: render one BashOperator's command of apogee_daily
with stand-in task context and print it. Runs in a child process so that the
environment (ALMANAC_BIN, AR_REPO, ...) seen by ar_common at import time can
differ per test.

    python _render_task.py TASK_ID OUTDIR MJD [key=value params ...]
"""

import os
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dags"))
warnings.simplefilter("ignore")


class FakeTI:
    """Answers the only xcom the bash templates pull: <tele>.resolve_night."""

    def __init__(self, outdir, mjd):
        self.outdir, self.mjd = outdir, mjd

    def xcom_pull(self, task_ids=None, **_):
        import ar_common as C
        tele = task_ids.split(".")[0]
        return C.night_paths(tele, self.mjd, self.outdir)


def params_of(dag):
    p = dag.params
    return dict(p.dump()) if hasattr(p, "dump") else dict(p)


def render(task, outdir, mjd, overrides=None):
    import jinja2
    params = params_of(task.dag)
    params.update(overrides or {})
    env = jinja2.Environment(undefined=jinja2.StrictUndefined)
    return env.from_string(task.bash_command).render(
        params=params, ti=FakeTI(outdir, mjd))


def main():
    task_id, outdir, mjd = sys.argv[1], sys.argv[2], int(sys.argv[3])
    overrides = dict(a.split("=", 1) for a in sys.argv[4:])
    import apogee_daily
    task = apogee_daily.dag.get_task(task_id)
    sys.stdout.write(render(task, outdir, mjd, overrides))


if __name__ == "__main__":
    main()
