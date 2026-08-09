import { apiGet } from "./client";

import type {
  AuditSummary,
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
