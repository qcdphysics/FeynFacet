#!/bin/bash
# Build and publish the generated pages on a dedicated GitHub Pages branch.
# The source code remains on its existing branch. Configure GitHub Pages to
# serve gh-pages (or the dedicated branch named by SITE_BRANCH), at its root.
#   SITE_REPO=git@github.com:qcdphysics/FeynFacet.git bash site/publish.sh
set -euo pipefail
# The website stays private until the paper is posted: no publication without the arXiv link.
if [ -z "${PAPER_URL:-}" ]; then
  echo "site/publish.sh: the site is published only after the paper is posted; set PAPER_URL to its arXiv link" >&2
  exit 2
fi
if [ -z "${SITE_REPO:-}" ]; then
  echo "site/publish.sh: set SITE_REPO to the Git repository that hosts GitHub Pages" >&2
  exit 2
fi
site_branch=${SITE_BRANCH:-gh-pages}
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
git check-ref-format --branch "$site_branch" >/dev/null
source_branch=$(git -C "$project_root" symbolic-ref --quiet --short HEAD || true)
case "$site_branch" in
  main|master)
    echo "site/publish.sh: refusing to publish generated pages over $site_branch" >&2
    exit 2 ;;
esac
if [ "$site_branch" = "$source_branch" ]; then
  echo "site/publish.sh: the Pages branch must differ from the current source branch" >&2
  exit 2
fi

# Record the exact remote head before building. An explicit lease also protects
# a new branch: an empty expected head means that the branch must still be absent.
remote_refs=$(git ls-remote --symref "$SITE_REPO" HEAD "refs/heads/$site_branch")
default_ref=$(printf '%s\n' "$remote_refs" | awk '$1 == "ref:" && $3 == "HEAD" {print $2}')
if [ "$default_ref" = "refs/heads/$site_branch" ]; then
  echo "site/publish.sh: refusing to replace the remote default branch" >&2
  exit 2
fi
expected_head=$(printf '%s\n' "$remote_refs" | awk -v ref="refs/heads/$site_branch" '$2 == ref && $1 != "ref:" {print $1}')
source_commit=$(git -C "$project_root" rev-parse --short HEAD)
author_name=$(git -C "$project_root" config user.name)
author_email=$(git -C "$project_root" config user.email)
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
git -C "$work" init -q -b "$site_branch"
if [ -n "$expected_head" ]; then
  git -C "$work" fetch -q --depth=1 "$SITE_REPO" "refs/heads/$site_branch"
  if [ "$(git -C "$work" rev-parse FETCH_HEAD)" != "$expected_head" ]; then
    echo "site/publish.sh: the remote Pages branch changed before the build; retry" >&2
    exit 2
  fi
  for source_path in FeynFacet/FeynFacet.wl Processes Start.wl; do
    if git -C "$work" cat-file -e "FETCH_HEAD:$source_path" 2>/dev/null; then
      echo "site/publish.sh: refusing to replace a branch that contains the project source" >&2
      exit 2
    fi
  done
fi

cd "$project_root"
python3 site/build.py
cp -R site/public/. "$work/"
cat > "$work/README.md" <<'TXT'
# FeynFacet website

The generated pages of the FeynFacet results website, published with GitHub
Pages. They are built from the source code by `site/build.py`; do not edit them
here.
TXT
touch "$work/.nojekyll"
git -C "$work" add -A
git -C "$work" -c user.name="$author_name" -c user.email="$author_email" \
  commit -q -m "Website built from the FeynFacet sources ($source_commit)"
git -C "$work" push -q \
  --force-with-lease="refs/heads/$site_branch:$expected_head" \
  "$SITE_REPO" "HEAD:refs/heads/$site_branch"
echo "published to $SITE_REPO ($site_branch)"
