# ExonDomainCompare frontend

The React frontend is normally installed and started from the repository root:

```bash
./scripts/setup_local.sh
./scripts/start_local.sh
```

Frontend-only development commands:

```bash
npm run dev
npm run lint
npm run parity
npm run build
```

The frontend expects the FastAPI backend at `http://127.0.0.1:8000` unless
`VITE_API_BASE` is configured.
