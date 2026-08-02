import { useState } from "react";

export default function CopyButton({ text }) {
  const [copied, setCopied] = useState(false);
  async function copy() {
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      // fallback for non-secure contexts
      const ta = document.createElement("textarea");
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      try { document.execCommand("copy"); } catch { /* ignore */ }
      document.body.removeChild(ta);
    }
    setCopied(true);
    setTimeout(() => setCopied(false), 1600);
  }
  return (
    <button className={`copy-btn${copied ? " copied" : ""}`} onClick={copy} type="button">
      {copied ? "Copied ✓" : "Copy"}
    </button>
  );
}
