// Fail-closed handoff adapter: never executes inherited foreign campaign commands.
// Current release also remains open at the explicit source-merge approval boundary.
const ledger = process.argv.at(-1) || '(no ledger)';
console.error(`UNMET: ${ledger}: historical lanes are outside this task; current release awaits source merge, successor artifact and release approval.`);
process.exit(1);
