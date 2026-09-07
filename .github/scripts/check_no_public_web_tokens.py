#!/usr/bin/env python3
"""Fail if any tracked file hard-codes a personal public-web URL.

Flatiron serves each user's public web directory at

    https://users.flatironinstitute.org/~<user>/<TOKEN>/...

where <TOKEN> is a long random path segment that acts as an unauthenticated
bearer credential: anyone who learns it can read that whole directory. These
repositories are public, so committing such a URL publishes the credential.

This check flags a URL whose segment after ~<user>/ is a long literal string.
It deliberately does NOT flag the URL built at runtime from an environment
variable (e.g. `.../~asaydjari/$(ENV["SLACK_TOKEN"])/...`), which is the
correct pattern.

Usage:
    python3 .github/scripts/check_no_public_web_tokens.py [paths...]
With no arguments it checks every tracked file in the repository.
Also usable as a pre-commit hook.
"""
import re
import subprocess
import sys

# Host + ~user + a long literal path segment == a hard-coded access token.
PATTERN = re.compile(r'users\.flatironinstitute\.org/~[A-Za-z0-9._-]+/[A-Za-z0-9_-]{16,}')

SKIP_SUFFIXES = ('.png', '.jpg', '.jpeg', '.gif', '.pdf', '.h5', '.hdf5',
                 '.fits', '.jld2', '.gz', '.zip', '.ico')


def tracked_files():
    out = subprocess.run(['git', 'ls-files', '-z'],
                         capture_output=True, text=True, check=True).stdout
    return [f for f in out.split('\0') if f]


def main(argv):
    paths = argv[1:] or tracked_files()
    # Never report this file's own explanatory text as a finding.
    self_name = 'check_no_public_web_tokens.py'
    findings = []
    for path in paths:
        if path.lower().endswith(SKIP_SUFFIXES) or path.endswith(self_name):
            continue
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
                for lineno, line in enumerate(fh, 1):
                    if PATTERN.search(line):
                        findings.append((path, lineno))
        except (IsADirectoryError, FileNotFoundError, PermissionError):
            continue

    if not findings:
        print(f'OK: no hard-coded public-web URLs in {len(paths)} tracked files.')
        return 0

    print('FAIL: hard-coded personal public-web URL(s) found.')
    print('These path segments are unauthenticated access tokens for a')
    print('personal web directory, and this repository is public.\n')
    for path, lineno in findings:
        # Print the location only -- never the matched text, which is the secret.
        print(f'  {path}:{lineno}')
    print('\nFix: reference the figure or directory by its internal path, e.g.')
    print('  /mnt/ceph/users/<user>/working/<date>/plots/<name>/')
    print('To show a figure in a PR, attach it instead of linking it:')
    print('  gh pr comment <n> -R <owner>/<repo> --attach fig.png --body "caption"')
    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
