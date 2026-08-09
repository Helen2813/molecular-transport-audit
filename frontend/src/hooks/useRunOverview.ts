import {
  useEffect,
  useState,
} from "react";

import {
  getRun,
  getRunSummary,
} from "../api/runs";

import type {
  AuditSummary,
  RunDetail,
} from "../types/audit";

interface RunOverviewState {
  run: RunDetail | null;
  summary: AuditSummary | null;

  loading: boolean;
  error: string | null;
}

export function useRunOverview(
  runId: string,
): RunOverviewState {
  const [
    state,
    setState,
  ] = useState<RunOverviewState>({
    run: null,
    summary: null,
    loading: true,
    error: null,
  });

  useEffect(() => {
    const controller =
      new AbortController();

    async function load() {
      setState((current) => ({
        ...current,
        loading: true,
        error: null,
      }));

      try {
        const [
          run,
          summary,
        ] = await Promise.all([
          getRun(
            runId,
            controller.signal,
          ),

          getRunSummary(
            runId,
            controller.signal,
          ),
        ]);

        setState({
          run,
          summary,
          loading: false,
          error: null,
        });
      } catch (error) {
        if (
          error instanceof DOMException &&
          error.name === "AbortError"
        ) {
          return;
        }

        setState({
          run: null,
          summary: null,
          loading: false,
          error:
            error instanceof Error
              ? error.message
              : "Could not load audit.",
        });
      }
    }

    void load();

    return () => {
      controller.abort();
    };
  }, [runId]);

  return state;
}
