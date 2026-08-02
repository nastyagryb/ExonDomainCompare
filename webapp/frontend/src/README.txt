Frontend v7: FastAPI-connected version

Copy these files to:
  bachelorthesis/webapp/frontend/src/App.jsx
  bachelorthesis/webapp/frontend/src/App.css
  bachelorthesis/webapp/frontend/src/index.css

Start backend first:
  cd webapp/backend
  python -m uvicorn main:app --reload --port 8000

Start frontend:
  cd webapp/frontend
  npm run dev

Open:
  http://localhost:5173

The frontend reads data from:
  http://localhost:8000/api/dashboard
  http://localhost:8000/api/species
  http://localhost:8000/api/transcripts
  http://localhost:8000/api/workflow
  http://localhost:8000/api/downloads

Optional: create .env in webapp/frontend with:
  VITE_API_BASE=http://localhost:8000
