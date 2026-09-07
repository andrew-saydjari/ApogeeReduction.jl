# Contributing

## Sharing figures and results in issues and PRs

**This repository is public.** Anything written into a PR body, a comment, or a
tracked file is world-readable.

### Never link a personal public-web directory

Flatiron serves each user's public web directory at

```
https://users.flatironinstitute.org/~<user>/<TOKEN>/...
```

The `<TOKEN>` segment is not a folder name. It is an unauthenticated bearer
credential: anyone who reads it gets access to that user's entire public web
directory, with no login. Pasting such a URL into a PR, a comment, or a
committed file publishes that credential.

This applies to the rendered URL only. Building the URL at runtime from an
environment variable is fine and is the established pattern:

```julia
"https://users.flatironinstitute.org/~asaydjari/$(ENV["SLACK_TOKEN"])/sdsswork/"
```

CI enforces this. `.github/workflows/secret-scan.yml` runs
`.github/scripts/check_no_public_web_tokens.py` on every push and PR and fails
if a tracked file contains a hard-coded one. You can run it locally at any time:

```bash
python3 .github/scripts/check_no_public_web_tokens.py
```

### Do this instead

**To show a figure in a PR or issue: attach it.** Attachments are uploaded to
GitHub's own `user-attachments` store. They are visible in the thread, they do
not require committing binaries to the repository, and they leak nothing.

```bash
gh pr comment <number> -R <owner>/<repo> \
  --attach ./fig1.png#"alt text" \
  --body "Caption explaining what the figure shows."
```

Notes, learned the hard way:

- Requires **`gh` 2.99.0 or newer** (`--attach` did not exist before that). The
  cluster Lmod `gh` is older; a current static binary is installed at
  `~/bin/gh` on the Flatiron side.
- Use **`gh pr comment --attach`**, not `gh pr edit --attach`. `gh pr edit`
  issues a metadata GraphQL query that needs the `read:org` scope, which our
  tokens do not carry, so it fails. `gh pr comment` makes no such query.
- If the `--body` markdown already references the file (``![alt](./fig1.png)``),
  that reference is rewritten in place to point at the uploaded asset, so you
  control where each image appears and can interleave captions. Any attached
  file the body does not reference is appended at the end.
- Up to 50 files per invocation.
- `gh issue comment` takes the same `--attach` flag.

**To point at a directory of outputs, or at a long write-up you do not want to
publish: use the internal path.**

```
/mnt/ceph/users/<user>/working/<YYYY_MM_DD>/plots/<name>/
```

That path is not sensitive. It is only reachable by someone who already has
cluster access, which is exactly the audience for a raw output directory.

**Do not commit figures into the repository** to work around this. Attach them.

### If a token URL has already been committed or posted

1. Edit the PR body, comment, or file to remove it, and open the fix as a PR.
2. Editing a PR body or comment is **not sufficient on its own.** GitHub keeps
   the previous revision, and anyone with read access to a public repo can read
   it through the "edited" dropdown. The old revision must also be deleted:
   open the comment, click **edited**, pick the revision, then
   **Options → Delete revision from history**. There is no API for this; it is
   web-UI only, and only the author or someone with write access can do it.
3. A URL that reached a commit stays reachable at that commit even after a
   follow-up fix. Removing it fully requires rewriting or deleting that branch.
4. Actions run logs cannot be edited; the only remedy is deleting the run.

If any of those cannot be completed, treat the token as compromised and rotate
it rather than assuming the edit was enough.
