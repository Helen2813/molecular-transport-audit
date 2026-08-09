import {
  AlertTriangle,
  RefreshCw,
} from "lucide-react";

import OverviewSkeleton from "../components/overview/OverviewSkeleton";
import OverviewStats from "../components/overview/OverviewStats";
import RunHeader from "../components/overview/RunHeader";

import { useRunOverview } from "../hooks/useRunOverview";

const RUN_ID =
  "unified_audit_validation";

export default function OverviewPage() {
  const {
    run,
    summary,
    loading,
    error,
  } = useRunOverview(
    RUN_ID,
  );

  return (
    <main className="min-h-screen pl-[248px] pt-[72px]">
      <div className="min-h-[calc(100vh-72px)] px-7 pb-12 pt-6">
        <div className="mx-auto max-w-[1600px]">
          {loading && (
            <OverviewSkeleton />
          )}

          {!loading && error && (
            <div className="rounded-xl border border-red-400/15 bg-red-500/[0.035] p-5">
              <div className="flex items-start gap-3">
                <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-lg bg-red-500/10">
                  <AlertTriangle
                    size={18}
                    className="text-red-400"
                  />
                </div>

                <div>
                  <h2 className="text-sm font-semibold text-red-200">
                    Could not load audit
                  </h2>

                  <p className="mt-1 text-xs leading-5 text-red-200/60">
                    {error}
                  </p>

                  <button
                    type="button"
                    onClick={() =>
                      window.location.reload()
                    }
                    className="mt-4 flex items-center gap-2 rounded-lg border border-red-400/15 bg-red-400/[0.06] px-3 py-2 text-xs font-medium text-red-200 transition hover:bg-red-400/10"
                  >
                    <RefreshCw
                      size={14}
                    />

                    Retry
                  </button>
                </div>
              </div>
            </div>
          )}

          {!loading &&
            !error &&
            run &&
            summary && (
              <>
                <RunHeader
                  externalCohort={
                    summary.run
                      .external_cohort
                  }
                />

                <OverviewStats
                  run={run}
                  summary={summary}
                />

                <div className="mt-4 min-h-[520px] rounded-[14px] border border-dashed border-white/[0.04] bg-white/[0.008]">
                  <div className="flex h-24 items-center justify-center text-xs text-slate-600">
                    Module Results will go here
                  </div>
                </div>
              </>
            )}
        </div>
      </div>
    </main>
  );
}
