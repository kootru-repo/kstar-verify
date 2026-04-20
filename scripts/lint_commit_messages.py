"""Lint commit messages in a range for process-narrative vocabulary.

Usage:
  python3 scripts/lint_commit_messages.py <base-sha> <head-sha>

Exit code 0 if every commit in (base..head] has a clean subject and
body; 1 if any pattern matches.  Intended for PR-time CI so
breadcrumb-trail commit messages never land on main.

The pattern library is kept in one place and documented inline.  If a
match is legitimate (e.g., a paper-source edit that genuinely refers
to an external identifier), amend the commit to describe the change
in content terms instead.
"""
from __future__ import annotations

import re
import subprocess
import sys

# Each tuple: (label, regex, hint shown on failure).
PATTERNS = [
    # Cleanup / process narrative ----------------------------------------
    ("cleanup-verb-scrub",
     re.compile(r"\bscrub\w*", re.I),
     "describe the content change, not the cleanup verb"),
    ("cleanup-verb-strip",
     re.compile(r"\bstrip\w*", re.I),
     "describe what was replaced, not the act of removing"),
    ("cleanup-verb-tidy",
     re.compile(r"\btidy\b", re.I),
     "'tidy' signals a post-hoc polish pass; describe the content"),
    ("cleanup-verb-condense",
     re.compile(r"\bcondens\w*", re.I),
     "describe the content change, not the compression"),
    ("cleanup-verb-polish",
     re.compile(r"\bpolish\w*", re.I),
     "describe the content change, not the act of polishing"),
    ("cleanup-verb-generic",
     re.compile(r"\bmake\s+(?:repo|it|this)\s+(?:generic|agnostic)\b", re.I),
     "renaming narratives should be avoided in commit subjects"),
    ("cleanup-agnostic",
     re.compile(r"\b(?:submission|journal)[- ]agnostic\b", re.I),
     "'agnostic' is a rebranding narrative"),
    ("cleanup-refresh",
     re.compile(r"\brefresh\w*", re.I),
     "name the file/content being updated rather than 'refresh'"),
    # Iteration / pass markers ------------------------------------------
    ("pass-counter-nth",
     re.compile(r"\b(?:first|second|third|fourth|fifth)\s+pass\b", re.I),
     "'Nth pass' reveals multiple cleanup rounds"),
    ("pass-counter-remaining",
     re.compile(r"\bremaining\s+\w+\s+(?:counts?|refs?|references?|mentions?)\b", re.I),
     "'remaining X' implies prior incomplete pass"),
    # Review / referee vocabulary ---------------------------------------
    ("review-referee",
     re.compile(r"\brefere?e\b", re.I),
     "avoid naming the reviewer audience in commit messages"),
    ("review-bulletproof",
     re.compile(r"\bbulletproof\w*", re.I),
     "describe the technical change, not the review defence"),
    ("review-posturing",
     re.compile(r"\bposturing\b", re.I),
     "review-adversarial vocabulary"),
    ("review-defensive",
     re.compile(r"\bdefensive\s+(?:framing|signal|commentary)\b", re.I),
     "review-adversarial vocabulary"),
    # External model attributions ---------------------------------------
    ("llm-claude",
     re.compile(r"\b(?:claude|anthropic|opus)\b", re.I),
     "external model attribution"),
    ("llm-gpt",
     re.compile(r"\bGPT[- ]?[0-9]+(?:\.[0-9]+)?\b", re.I),
     "external model attribution"),
    ("llm-openai",
     re.compile(r"\b(?:ChatGPT|OpenAI)\b", re.I),
     "external model attribution"),
    ("llm-ai-use",
     re.compile(r"\bAI[- ](?:use|assisted|generated|review)\b", re.I),
     "AI-use narrative"),
    ("llm-large-model",
     re.compile(r"\b(?:large language model|LLM)\b"),
     "LLM attribution"),
    # Private tooling ---------------------------------------------------
    ("private-theorem-cert",
     re.compile(r"\btheorem[- ]cert\w*", re.I),
     "private tooling name"),
    ("private-coherence-engine",
     re.compile(r"\bcoherence\s+engine\b", re.I),
     "private tooling name"),
    ("private-claim-engine",
     re.compile(r"\bclaim\s+engine\b", re.I),
     "private tooling name"),
    ("private-predicate-engine",
     re.compile(r"\bpredicate\s+engine\b", re.I),
     "private tooling name"),
    ("private-paste-stub",
     re.compile(r"\bpaste[- ]stubs?\b", re.I),
     "private tooling concept"),
    ("private-triple-mirror",
     re.compile(r"\bthree[- ]location\s+sync\b", re.I),
     "private workflow concept"),
    # Journal specificity ----------------------------------------------
    ("journal-prxq",
     re.compile(r"\bPRX[- ]?Q\b|PRX Quantum|PRXQ|prxq", re.I),
     "journal-specific reference"),
    ("journal-jgp",
     re.compile(r"\bJGP\b"),
     "journal-specific reference"),
    ("journal-craton-theorem",
     re.compile(r"\bCraton\s+theorem\b", re.I),
     "name a theorem by its content, not by the author"),
    # Local workspace paths --------------------------------------------
    ("path-local-win",
     re.compile(r"C:\\\\(?:Users|Projects)"),
     "local filesystem path"),
    ("path-local-unix",
     re.compile(r"/home/[a-z]|/Users/[a-zA-Z]"),
     "local filesystem path"),
    ("path-private-workspace",
     re.compile(r"\bSpectral-Geometry\b|\bkstar-certification\b"),
     "private workspace directory name"),
]


def _git_log(base: str, head: str) -> list[tuple[str, str, str]]:
    """Return [(sha, subject, body), ...] for commits in base..head."""
    # Use git's own format escapes so no NUL byte is embedded in argv
    # (Windows subprocess rejects embedded nulls in argument strings).
    fmt = "%H%x1f%s%x1f%b%x00"
    out = subprocess.check_output(
        ["git", "log", f"{base}..{head}", "--format=" + fmt],
    ).decode("utf-8", errors="replace")
    sep = "\x1f"
    term = "\x00"
    commits: list[tuple[str, str, str]] = []
    for block in out.split(term):
        if not block.strip():
            continue
        parts = block.split(sep, 2)
        if len(parts) < 3:
            continue
        sha, subject, body = parts
        commits.append((sha.strip(), subject, body))
    return commits


def _scan(subject: str, body: str) -> list[tuple[str, str, str, str]]:
    hits: list[tuple[str, str, str, str]] = []
    for label, regex, hint in PATTERNS:
        for m in regex.finditer(subject):
            hits.append((label, m.group(0), "subject", hint))
        for m in regex.finditer(body):
            hits.append((label, m.group(0), "body", hint))
    return hits


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: lint_commit_messages.py <base-sha> <head-sha>",
              file=sys.stderr)
        return 2
    base, head = argv[1], argv[2]
    commits = _git_log(base, head)
    any_fail = False
    for sha, subject, body in commits:
        hits = _scan(subject, body)
        if not hits:
            continue
        any_fail = True
        print(f"\n::error::commit {sha[:8]} has disallowed terms in its message")
        print(f"  subject: {subject}")
        for label, matched, loc, hint in hits:
            print(f"    - [{label}] {matched!r} (in {loc}): {hint}")
    if any_fail:
        print(
            "\nOne or more commits contain process-narrative vocabulary.\n"
            "Amend the offending commits (e.g. `git rebase -i` + reword) "
            "and push again.  Describe the content change, not the process."
        )
        return 1
    print(f"OK: {len(commits)} commit(s) scanned, no disallowed terms.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
