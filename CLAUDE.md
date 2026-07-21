# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@AGENTS.md

Everything below is specific to Claude Code's own surfaces (conversation lifecycle,
scheduling tools, `~/.claude/` paths). It has no Codex equivalent and is deliberately
not in `AGENTS.md`.

## Critical: do not archive this project (added 2026-05-18)

Never archive, close, end, hide, delete, or mark this project/conversation
as complete unless the user types this exact phrase:

> **`Archive this project now.`**

Soft-sounding user phrases that do NOT mean archive: "wrap up", "summarize",
"finalize", "done", "pause", "stop". When you see any of those, write a
summary only — do not archive anything, do not end the conversation, do
not call any tool that ends the session.

This applies regardless of session lifecycle events (compaction, session
end, context cleanup). "Tidying", "housekeeping", and "cleanup" are not
autonomous actions in this project.

Concrete behaviours:

- **`ScheduleWakeup` is forbidden in this project.** It announces "Nothing
  more to do this turn" which the UI reads as dormancy and archives the
  conversation. Use foreground polling, `Monitor`, or a background `Bash`
  with `run_in_background` followed by an active wait instead.
- **Don't end turns with no tool calls when work is in flight.** Even a
  brief polling tool call beats an empty turn.
- **No autonomous file moves into `archive/`, `old/`, `deprecated/`, or
  similar directories** without an explicit, in-chat instruction from the
  user naming the specific files. (Investigation 2026-05-18 found no such
  mover wired up; this rule is here defensively for future agents.)

## Host-specific notes (optional)

On the primary author's machine, two persistent memory notes live under `~/.claude/projects/-Users-z3437171-Dropbox-Github-Local-pigauto/memory/`: `user_profile.md` (author priorities) and `project_bace.md` (BACE internals). They are **not** portable — ignore this section if the path does not exist on the current host.
