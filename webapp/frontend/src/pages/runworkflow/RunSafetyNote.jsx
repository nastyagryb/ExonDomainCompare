export default function RunSafetyNote({ compact = false }) {
  return (
    <div className={`safety-note${compact ? " compact" : ""}`}>
      <span className="safety-icon" aria-hidden>🔒</span>
      <div className="safety-text">
        <strong>Cluster annotation runs in your local terminal.</strong>
        <p>
          The webapp does <b>not</b> collect LRZ passwords, 2FA codes or SSH credentials,
          and it never submits jobs to SLURM automatically. It only creates local run
          folders and reads <code>status.json</code>. Cluster steps are executed by you via
          the copy-paste commands below; login and 2FA happen on your own machine.
        </p>
      </div>
    </div>
  );
}
