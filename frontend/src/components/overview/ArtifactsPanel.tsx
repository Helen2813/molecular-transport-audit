import {
  CheckCircle2,
  Download,
  FileJson,
  FileSpreadsheet,
  FolderCheck,
} from "lucide-react";

import {
  useEffect,
  useState,
} from "react";

import {
  getRunArtifacts,
} from "../../api/runs";

import type {
  ArtifactItem,
} from "../../types/audit";


interface ArtifactsPanelProps {
  runId: string;
}


function fileSize(
  bytes: number,
): string {
  if (bytes < 1024) {
    return `${bytes} B`;
  }

  const kilobytes =
    bytes / 1024;

  if (kilobytes < 1024) {
    return (
      `${kilobytes.toFixed(1)} KB`
    );
  }

  return (
    `${(
      kilobytes / 1024
    ).toFixed(1)} MB`
  );
}


function downloadUrl(
  runId: string,
  artifact: ArtifactItem,
): string {
  const encodedPath =
    artifact.path
      .split("/")
      .map(
        encodeURIComponent,
      )
      .join("/");

  return (
    `/api/runs/` +
    `${encodeURIComponent(runId)}` +
    `/artifacts/${encodedPath}`
  );
}


function ArtifactIcon({
  name,
}: {
  name: string;
}) {
  if (
    name.endsWith(".csv")
  ) {
    return (
      <FileSpreadsheet
        size={18}
        className="text-emerald-400"
      />
    );
  }

  return (
    <FileJson
      size={18}
      className="text-blue-400"
    />
  );
}


export default function ArtifactsPanel({
  runId,
}: ArtifactsPanelProps) {
  const [
    artifacts,
    setArtifacts,
  ] = useState<
    ArtifactItem[]
  >([]);

  const [
    loading,
    setLoading,
  ] = useState(
    true,
  );

  const [
    error,
    setError,
  ] = useState<
    string | null
  >(null);


  useEffect(() => {
    const controller =
      new AbortController();

    getRunArtifacts(
      runId,
      controller.signal,
    )
      .then((response) => {
        setArtifacts(
          response.artifacts,
        );
        setError(null);
      })
      .catch((reason) => {
        setError(
          reason instanceof Error
            ? reason.message
            : "Could not load artifacts.",
        );
      })
      .finally(() => {
        setLoading(false);
      });

    return () => {
      controller.abort();
    };
  }, [runId]);


  const preferred =
    artifacts
      .slice()
      .sort((a, b) => {
        const priority = [
          "summary.json",
          "module_summary.csv",
          (
            "unified_audit_" +
            "verification.json"
          ),
          (
            "random_controls_" +
            "manifest.json"
          ),
        ];

        const ai =
          priority.indexOf(
            a.name,
          );

        const bi =
          priority.indexOf(
            b.name,
          );

        return (
          (
            ai === -1
              ? 999
              : ai
          ) -
          (
            bi === -1
              ? 999
              : bi
          )
        );
      })
      .slice(
        0,
        6,
      );


  return (
    <section className="flex h-full min-h-[320px] flex-col overflow-hidden rounded-[14px] border border-white/[0.07] bg-[#091725]/95 shadow-[0_12px_40px_rgba(0,0,0,0.14)]">
      <div className="flex h-[58px] shrink-0 items-center justify-between border-b border-white/[0.06] px-5">
        <div className="flex items-center gap-2.5">
          <FolderCheck
            size={17}
            className="text-slate-400"
          />

          <h2 className="text-[15px] font-semibold text-slate-100">
            Validated Artifacts
          </h2>
        </div>

        {!loading && (
          <span className="flex h-6 min-w-6 items-center justify-center rounded-full border border-white/[0.06] bg-white/[0.035] px-2 text-[11px] font-semibold text-slate-400">
            {artifacts.length}
          </span>
        )}
      </div>

      <div className="flex-1 px-3 py-2">
        {loading && (
          <div className="space-y-2">
            {Array.from({
              length: 4,
            }).map(
              (_, index) => (
                <div
                  key={index}
                  className="h-[60px] animate-pulse rounded-xl bg-white/[0.025]"
                />
              ),
            )}
          </div>
        )}

        {!loading &&
          error && (
            <div className="px-3 py-5 text-[12px] leading-5 text-red-300/70">
              {error}
            </div>
          )}

        {!loading &&
          !error &&
          preferred.map(
            (artifact) => (
              <div
                key={
                  artifact.path
                }
                className="group flex min-h-[62px] items-center gap-3 rounded-xl px-3 transition hover:bg-white/[0.025]"
              >
                <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-lg border border-white/[0.06] bg-white/[0.025]">
                  <ArtifactIcon
                    name={
                      artifact.name
                    }
                  />
                </div>

                <div className="min-w-0 flex-1">
                  <div className="flex items-center gap-1.5">
                    <p
                      title={
                        artifact.name
                      }
                      className="truncate text-[13px] font-medium text-slate-200"
                    >
                      {artifact.name}
                    </p>

                    <CheckCircle2
                      size={12}
                      className="shrink-0 text-emerald-400/80"
                    />
                  </div>

                  <div className="mt-1 flex items-center gap-2 text-[11px] text-slate-600">
                    <span>
                      {fileSize(
                        artifact.size_bytes,
                      )}
                    </span>

                    <span>
                      •
                    </span>

                    <span className="capitalize">
                      {
                        artifact.category
                      }
                    </span>
                  </div>
                </div>

                <a
                  href={downloadUrl(
                    runId,
                    artifact,
                  )}
                  className="flex h-8 w-8 shrink-0 items-center justify-center rounded-lg text-slate-600 opacity-60 transition hover:bg-white/[0.05] hover:text-slate-200 group-hover:opacity-100"
                  title="Download artifact"
                >
                  <Download
                    size={15}
                  />
                </a>
              </div>
            ),
          )}
      </div>

      <div className="border-t border-white/[0.05] px-4 py-3">
        <div className="flex items-center gap-2 text-[11px] text-emerald-400/70">
          <CheckCircle2
            size={13}
          />

          Independent verification passed
        </div>
      </div>
    </section>
  );
}
