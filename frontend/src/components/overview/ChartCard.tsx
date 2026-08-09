import {
  Info,
  type LucideIcon,
} from "lucide-react";

import type {
  ReactNode,
} from "react";


interface ChartCardProps {
  title: string;
  subtitle?: string;
  icon: LucideIcon;
  children: ReactNode;
  action?: ReactNode;
}


export default function ChartCard({
  title,
  subtitle,
  icon: Icon,
  children,
  action,
}: ChartCardProps) {
  return (
    <section className="group relative min-w-0 overflow-hidden rounded-[16px] border border-white/[0.07] bg-[#091725]/95 shadow-[0_18px_55px_rgba(0,0,0,0.16)]">
      <div className="pointer-events-none absolute inset-x-0 top-0 h-px bg-gradient-to-r from-transparent via-cyan-300/[0.10] to-transparent" />

      <div className="flex min-h-[66px] items-center justify-between border-b border-white/[0.06] px-5">
        <div className="flex min-w-0 items-center gap-3">
          <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-[10px] border border-white/[0.06] bg-white/[0.025]">
            <Icon
              size={17}
              strokeWidth={1.8}
              className="text-slate-300"
            />
          </div>

          <div className="min-w-0">
            <div className="flex items-center gap-2">
              <h3 className="truncate text-[15px] font-semibold tracking-[-0.015em] text-slate-100">
                {title}
              </h3>

              <Info
                size={13}
                className="shrink-0 text-slate-600"
              />
            </div>

            {subtitle && (
              <p className="mt-1 truncate text-[12px] text-slate-500">
                {subtitle}
              </p>
            )}
          </div>
        </div>

        {action}
      </div>

      <div className="px-4 pb-4 pt-5">
        {children}
      </div>
    </section>
  );
}
