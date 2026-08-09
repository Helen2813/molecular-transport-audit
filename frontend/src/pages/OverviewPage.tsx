import {
  AlertTriangle,
  LoaderCircle,
  RefreshCw,
} from "lucide-react";

import ModuleResultsTable from "../components/overview/ModuleResultsTable";
import OverviewCharts from "../components/overview/OverviewCharts";
import OverviewSkeleton from "../components/overview/OverviewSkeleton";
import OverviewStats from "../components/overview/OverviewStats";
import RunHeader from "../components/overview/RunHeader";
import RunPipeline from "../components/overview/RunPipeline";
import ArtifactsPanel from "../components/overview/ArtifactsPanel";

import {
  useRunOverview,
} from "../hooks/useRunOverview";


const DEFAULT_RUN_ID =
  "unified_audit_validation";


function resolveRunId(): string {
  const params =
    new URLSearchParams(
      window.location.search,
    );

  return (
    params.get("run") ??
    DEFAULT_RUN_ID
  );
}


export default function OverviewPage() {
  const runId =
    resolveRunId();

  const {
    run,
    summary,
    loading,
    error,
  } = useRunOverview(
    runId,
  );

  return (
    <main className="min-h-screen pl-[248px] pt-[72px]">
      <div className="min-h-[calc(100vh-72px)] px-7 pb-12 pt-6">
        <div className="mx-auto max-w-[1600px]">
          {loading && (
            <OverviewSkeleton />
          )}

          {!loading &&
            error && (
              <div className="rounded-xl border border-red-400/15 bg-red-500/[0.035] p-5">
                <div className="flex items-start gap-3">
                  <div className="flex h-10 w-10 shrink-0 items-center justify-center rounded-lg bg-red-500/10">
                    <AlertTriangle
                      size={19}
                      className="text-red-400"
                    />
                  </div>

                  <div>
                    <h2 className="text-[15px] font-semibold text-red-200">
                      Could not load audit
                    </h2>

                    <p className="mt-1 text-[13px] leading-5 text-red-200/60">
                      {error}
                    </p>

                    <button
                      type="button"
                      onClick={() =>
                        window.location.reload()
                      }
                      className="mt-4 flex items-center gap-2 rounded-lg border border-red-400/15 bg-red-400/[0.06] px-3 py-2 text-[12px] font-medium text-red-200 transition hover:bg-red-400/10"
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
            run && (
              <>
                <RunHeader
                  externalCohort={
                    summary?.run
                      .external_cohort ??
                    run.external_cohort ??
                    "External"
                  }
                />

                <OverviewStats
                  run={run}
                  summary={summary}
                />

                {!summary && (
                  <>
                    <div className="mt-4">
                      <RunPipeline
                        run={run}
                      />
                    </div>

                    {run.status !==
                      "failed" && (
                      <div className="mt-4 flex min-h-[150px] items-center justify-center rounded-[16px] border border-dashed border-white/[0.06] bg-white/[0.01]">
                        <div className="text-center">
                          <div className="mx-auto flex h-10 w-10 items-center justify-center rounded-full border border-cyan-400/10 bg-cyan-400/[0.05]">
                            <LoaderCircle
                              size={18}
                              className="animate-spin text-cyan-400"
                            />
                          </div>

                          <p className="mt-3 text-[14px] font-medium text-slate-300">
                            Audit results are being generated
                          </p>

                          <p className="mt-1 text-[12px] text-slate-600">
                            Results become available only after successful verification.
                          </p>
                        </div>
                      </div>
                    )}
                  </>
                )}

                {summary && (
                  <>
                    <div className="mt-4 grid grid-cols-1 gap-4 2xl:grid-cols-[minmax(0,1fr)_330px]">
                      <ModuleResultsTable
                        modules={summary.modules}
                      />

                      <ArtifactsPanel
                        runId={run.run_id}
                      />
                    </div>

                    <div className="mt-4">
                      <OverviewCharts
                        modules={
                          summary.modules
                        }
                      />
                    </div>

                    <div className="mt-4">
                      <RunPipeline
                        run={run}
                      />
                    </div>
                  </>
                )}
              </>
            )}
        </div>
      </div>
    </main>
  );
}
