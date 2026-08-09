import {
  useEffect,
  useState,
} from "react";

import {
  getRun,
  getRunSummary,
} from "../api/runs";

import {
  ApiError,
} from "../api/client";

import type {
  AuditSummary,
  RunDetail,
} from "../types/audit";


interface RunOverviewState {
  run: RunDetail | null;
  summary: AuditSummary | null;

  loading: boolean;
  refreshing: boolean;

  error: string | null;
}


const POLLING_INTERVAL_MS =
  2500;


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
    refreshing: false,
    error: null,
  });

  useEffect(() => {
    let cancelled = false;

    let timer:
      | number
      | undefined;

    let controller:
      | AbortController
      | null = null;


    async function load(
      initial: boolean,
    ) {
      controller?.abort();

      controller =
        new AbortController();

      if (!initial) {
        setState(
          (current) => ({
            ...current,
            refreshing: true,
          }),
        );
      }

      try {
        const run =
          await getRun(
            runId,
            controller.signal,
          );

        let summary:
          | AuditSummary
          | null = null;

        if (
          run.status ===
            "completed" &&
          run.validated
        ) {
          try {
            summary =
              await getRunSummary(
                runId,
                controller.signal,
              );
          } catch (error) {
            if (
              error instanceof
                ApiError &&
              error.status === 409
            ) {
              summary = null;
            } else {
              throw error;
            }
          }
        }

        if (cancelled) {
          return;
        }

        setState({
          run,
          summary,
          loading: false,
          refreshing: false,
          error: null,
        });

        if (
          run.status ===
            "queued" ||
          run.status ===
            "running"
        ) {
          timer =
            window.setTimeout(
              () => {
                void load(
                  false,
                );
              },
              POLLING_INTERVAL_MS,
            );
        }
      } catch (error) {
        if (
          cancelled ||
          (
            error instanceof
              DOMException &&
            error.name ===
              "AbortError"
          )
        ) {
          return;
        }

        setState(
          (current) => ({
            ...current,
            loading: false,
            refreshing: false,
            error:
              error instanceof
                Error
                ? error.message
                : (
                  "Could not " +
                  "load audit."
                ),
          }),
        );
      }
    }


    void load(
      true,
    );


    return () => {
      cancelled = true;

      controller?.abort();

      if (
        timer !== undefined
      ) {
        window.clearTimeout(
          timer,
        );
      }
    };
  }, [runId]);


  return state;
}
