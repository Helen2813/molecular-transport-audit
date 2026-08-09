const API_BASE_URL =
  import.meta.env.VITE_API_BASE_URL ??
  "/api";


export class ApiError extends Error {
  readonly status: number;

  constructor(
    message: string,
    status: number,
  ) {
    super(message);

    this.name = "ApiError";
    this.status = status;
  }
}


async function apiRequest<T>(
  path: string,
  options: RequestInit,
): Promise<T> {
  const response = await fetch(
    `${API_BASE_URL}${path}`,
    {
      ...options,

      headers: {
        Accept: "application/json",
        ...(options.body
          ? {
              "Content-Type":
                "application/json",
            }
          : {}),
        ...options.headers,
      },
    },
  );

  if (!response.ok) {
    let message =
      `API request failed: ${response.status}`;

    try {
      const payload =
        await response.json();

      if (
        typeof payload.detail ===
        "string"
      ) {
        message =
          payload.detail;
      } else if (
        payload.detail &&
        typeof payload.detail ===
          "object" &&
        typeof payload.detail
          .message === "string"
      ) {
        message =
          payload.detail.message;
      }
    } catch {
      // Keep generic HTTP error.
    }

    throw new ApiError(
      message,
      response.status,
    );
  }

  return response.json() as Promise<T>;
}


export function apiGet<T>(
  path: string,
  signal?: AbortSignal,
): Promise<T> {
  return apiRequest<T>(
    path,
    {
      method: "GET",
      signal,
    },
  );
}


export function apiPost<T>(
  path: string,
  body: unknown,
  signal?: AbortSignal,
): Promise<T> {
  return apiRequest<T>(
    path,
    {
      method: "POST",
      body: JSON.stringify(
        body,
      ),
      signal,
    },
  );
}
