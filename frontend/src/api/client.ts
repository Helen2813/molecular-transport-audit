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

export async function apiGet<T>(
  path: string,
  signal?: AbortSignal,
): Promise<T> {
  const response = await fetch(
    `${API_BASE_URL}${path}`,
    {
      method: "GET",

      headers: {
        Accept: "application/json",
      },

      signal,
    },
  );

  if (!response.ok) {
    let message =
      `API request failed: ${response.status}`;

    try {
      const payload = await response.json();

      if (typeof payload.detail === "string") {
        message = payload.detail;
      }
    } catch {
      // Keep the generic HTTP error.
    }

    throw new ApiError(
      message,
      response.status,
    );
  }

  return response.json() as Promise<T>;
}
