import {
  apiGet,
  apiPost,
} from "./client";

import type {
  ArtifactListResponse,
  AuditSummary,
  RunCreateResponse,
  RunDetail,
} from "../types/audit";


export function getRun(
  runId: string,
  signal?: AbortSignal,
): Promise<RunDetail> {
  return apiGet<RunDetail>(
    `/runs/${encodeURIComponent(runId)}`,
    signal,
  );
}


export function getRunSummary(
  runId: string,
  signal?: AbortSignal,
): Promise<AuditSummary> {
  return apiGet<AuditSummary>(
    `/runs/${encodeURIComponent(runId)}/summary`,
    signal,
  );
}


export function getRunArtifacts(
  runId: string,
  signal?: AbortSignal,
): Promise<ArtifactListResponse> {
  return apiGet<ArtifactListResponse>(
    `/runs/${encodeURIComponent(runId)}/artifacts`,
    signal,
  );
}


export function createRun():
Promise<RunCreateResponse> {
  return apiPost<RunCreateResponse>(
    "/runs",
    {
      analysis:
        "gse239948_external_canine_audit",
    },
  );
}
