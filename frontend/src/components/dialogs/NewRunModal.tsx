import {
  AlertCircle,
  ArrowRight,
  CheckCircle2,
  Database,
  LoaderCircle,
  Network,
  ShieldCheck,
  X,
} from "lucide-react";

import {
  useEffect,
  useState,
} from "react";

import {
  createPortal,
} from "react-dom";

import {
  createRun,
} from "../../api/runs";


interface NewRunModalProps {
  open: boolean;
  onClose: () => void;
}


export default function NewRunModal({
  open,
  onClose,
}: NewRunModalProps) {
  const [
    submitting,
    setSubmitting,
  ] = useState(false);

  const [
    error,
    setError,
  ] = useState<string | null>(
    null,
  );


  useEffect(() => {
    console.log(
      "[NewRunModal] open changed:",
      open,
    );

    if (!open) {
      return;
    }

    const previousOverflow =
      document.body.style.overflow;

    document.body.style.overflow =
      "hidden";


    function handleKeyDown(
      event: KeyboardEvent,
    ) {
      if (
        event.key ===
          "Escape" &&
        !submitting
      ) {
        onClose();
      }
    }


    window.addEventListener(
      "keydown",
      handleKeyDown,
    );


    return () => {
      document.body.style.overflow =
        previousOverflow;

      window.removeEventListener(
        "keydown",
        handleKeyDown,
      );
    };
  }, [
    open,
    submitting,
    onClose,
  ]);


  if (!open) {
    return null;
  }


  console.log(
    "[NewRunModal] rendering PORTAL",
  );


  async function startRun() {
    console.log(
      "[NewRunModal] Start Audit clicked",
    );

    setSubmitting(true);
    setError(null);

    try {
      const created =
        await createRun();

      console.log(
        "[NewRunModal] run created:",
        created,
      );

      window.location.assign(
        `/?run=${encodeURIComponent(
          created.run_id,
        )}`,
      );
    } catch (reason) {
      console.error(
        "[NewRunModal] createRun failed:",
        reason,
      );

      setError(
        reason instanceof Error
          ? reason.message
          : "Could not start audit.",
      );

      setSubmitting(false);
    }
  }


  return createPortal(
    <div
      className="fixed inset-0 flex items-center justify-center p-6"
      style={{
        zIndex: 2147483647,
      }}
    >
      {/* Backdrop */}
      <div
        className="absolute inset-0 bg-[#020812]/80 backdrop-blur-[6px]"
        onClick={() => {
          if (!submitting) {
            onClose();
          }
        }}
      />

      {/* Decorative glow */}
      <div className="pointer-events-none absolute left-1/2 top-1/2 h-[520px] w-[720px] -translate-x-1/2 -translate-y-1/2 rounded-full bg-blue-500/[0.06] blur-[100px]" />

      {/* Modal */}
      <div
        className="
          relative
          z-10
          w-full
          max-w-[620px]
          overflow-hidden
          rounded-[22px]
          border
          border-white/[0.11]
          bg-[#091725]
          shadow-[0_35px_120px_rgba(0,0,0,0.72)]
        "
      >
        <div className="pointer-events-none absolute inset-x-0 top-0 h-px bg-gradient-to-r from-transparent via-cyan-300/40 to-transparent" />

        {/* Header */}
        <div className="flex items-start justify-between border-b border-white/[0.07] px-7 py-6">
          <div>
            <p className="mb-2 text-[11px] font-semibold uppercase tracking-[0.14em] text-cyan-400">
              Scientific Workflow
            </p>

            <h2 className="text-[22px] font-semibold tracking-[-0.03em] text-white">
              New Molecular Transport Audit
            </h2>

            <p className="mt-2 text-[13px] leading-5 text-slate-500">
              Launch the validated external
              canine representation workflow.
            </p>
          </div>

          <button
            type="button"
            disabled={submitting}
            onClick={onClose}
            className="flex h-10 w-10 shrink-0 items-center justify-center rounded-xl border border-transparent text-slate-500 transition hover:border-white/[0.06] hover:bg-white/[0.04] hover:text-white disabled:opacity-40"
          >
            <X size={19} />
          </button>
        </div>

        {/* Body */}
        <div className="space-y-4 px-7 py-6">
          <div className="rounded-[15px] border border-cyan-400/15 bg-gradient-to-br from-cyan-400/[0.07] via-blue-500/[0.035] to-transparent p-5">
            <div className="flex gap-4">
              <div className="flex h-11 w-11 shrink-0 items-center justify-center rounded-xl border border-cyan-400/15 bg-cyan-400/[0.08]">
                <Network
                  size={21}
                  className="text-cyan-300"
                />
              </div>

              <div>
                <div className="flex items-center gap-2">
                  <p className="text-[15px] font-semibold text-slate-100">
                    GSE239948 External Canine
                    Audit
                  </p>

                  <span className="rounded-md border border-emerald-400/15 bg-emerald-400/[0.06] px-2 py-0.5 text-[9px] font-semibold uppercase tracking-[0.08em] text-emerald-400">
                    Validated
                  </span>
                </div>

                <p className="mt-2 text-[12px] leading-5 text-slate-500">
                  Direct preservation,
                  permutation inference,
                  reliability, matched random
                  controls, summary assembly
                  and independent verification.
                </p>
              </div>
            </div>
          </div>

          {/* Cohorts */}
          <div className="grid grid-cols-2 gap-3">
            <div className="rounded-[13px] border border-white/[0.065] bg-white/[0.018] p-4">
              <div className="flex items-center gap-2 text-[12px] text-slate-500">
                <Database size={15} />
                Reference Cohort
              </div>

              <p className="mt-2.5 text-[15px] font-semibold tracking-[-0.015em] text-slate-200">
                GSE238110_DOG2
              </p>
            </div>

            <div className="rounded-[13px] border border-white/[0.065] bg-white/[0.018] p-4">
              <div className="flex items-center gap-2 text-[12px] text-slate-500">
                <Database size={15} />
                External Cohort
              </div>

              <p className="mt-2.5 text-[15px] font-semibold tracking-[-0.015em] text-slate-200">
                GSE239948
              </p>
            </div>
          </div>

          {/* Guard */}
          <div className="flex items-center justify-between rounded-[13px] border border-emerald-400/10 bg-emerald-400/[0.025] px-4 py-3.5">
            <div className="flex items-center gap-3">
              <ShieldCheck
                size={19}
                className="text-emerald-400"
              />

              <div>
                <p className="text-[13px] font-medium text-slate-200">
                  Outcome-blind execution
                </p>

                <p className="mt-1 text-[11px] text-slate-600">
                  No outcome data are loaded
                  by this workflow.
                </p>
              </div>
            </div>

            <CheckCircle2
              size={19}
              className="text-emerald-400"
            />
          </div>

          {/* Pipeline preview */}
          <div className="rounded-[13px] border border-white/[0.06] bg-black/[0.08] px-4 py-4">
            <p className="text-[11px] font-medium uppercase tracking-[0.09em] text-slate-600">
              Workflow
            </p>

            <div className="mt-3 flex items-center">
              {[
                "Preservation",
                "Inference",
                "Reliability",
                "Controls",
                "Verification",
              ].map(
                (label, index, values) => (
                  <div
                    key={label}
                    className="flex min-w-0 flex-1 items-center"
                  >
                    <div className="min-w-0">
                      <div className="h-2 w-2 rounded-full bg-cyan-400" />

                      <span className="mt-2 block truncate pr-2 text-[10px] text-slate-500">
                        {label}
                      </span>
                    </div>

                    {index <
                      values.length -
                        1 && (
                      <div className="mb-5 mr-2 h-px flex-1 bg-gradient-to-r from-cyan-400/40 to-white/[0.06]" />
                    )}
                  </div>
                ),
              )}
            </div>
          </div>

          {error && (
            <div className="flex gap-3 rounded-[13px] border border-red-400/15 bg-red-400/[0.045] px-4 py-3.5">
              <AlertCircle
                size={18}
                className="mt-0.5 shrink-0 text-red-400"
              />

              <div>
                <p className="text-[12px] font-semibold text-red-300">
                  Could not start audit
                </p>

                <p className="mt-1 text-[12px] leading-5 text-red-200/65">
                  {error}
                </p>
              </div>
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="flex items-center justify-between border-t border-white/[0.07] bg-black/[0.08] px-7 py-5">
          <p className="max-w-[285px] text-[11px] leading-5 text-slate-600">
            One scientific workflow can run
            at a time in the current backend.
          </p>

          <div className="flex items-center gap-3">
            <button
              type="button"
              disabled={submitting}
              onClick={onClose}
              className="h-11 rounded-xl border border-white/[0.07] px-5 text-[13px] font-medium text-slate-400 transition hover:bg-white/[0.035] hover:text-white disabled:opacity-40"
            >
              Cancel
            </button>

            <button
              type="button"
              disabled={submitting}
              onClick={() => {
                void startRun();
              }}
              className="flex h-11 min-w-[145px] items-center justify-center gap-2 rounded-xl border border-blue-300/15 bg-gradient-to-r from-[#2463f3] to-[#7651ef] px-5 text-[13px] font-semibold text-white shadow-[0_10px_30px_rgba(67,73,255,0.25)] transition hover:brightness-110 disabled:cursor-wait disabled:opacity-65"
            >
              {submitting ? (
                <>
                  <LoaderCircle
                    size={17}
                    className="animate-spin"
                  />
                  Starting...
                </>
              ) : (
                <>
                  Start Audit
                  <ArrowRight
                    size={16}
                  />
                </>
              )}
            </button>
          </div>
        </div>
      </div>
    </div>,
    document.body,
  );
}
