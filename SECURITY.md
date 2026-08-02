# Security

Do not report passwords, MFA responses, private keys, host-specific paths or
unredacted cluster output in a public issue.

Before sharing diagnostics, use:

```bash
.venv/bin/edc doctor --redact-paths
.venv/bin/edc cluster doctor --redact-paths
```

The application does not store LRZ passwords or MFA responses. Review every
remote installation plan before adding `--confirm`.

