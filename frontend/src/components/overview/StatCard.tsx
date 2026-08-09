import type {
  LucideIcon,
} from "lucide-react";

interface StatCardProps {
  label: string;
  value: string;
  description: string;

  icon: LucideIcon;

  iconClassName: string;
  iconBackgroundClassName: string;

  valueClassName?: string;
}

export default function StatCard({
  label,
  value,
  description,
  icon: Icon,
  iconClassName,
  iconBackgroundClassName,
  valueClassName = "text-white",
}: StatCardProps) {
  return (
    <div className="group relative min-w-0 overflow-hidden rounded-[13px] border border-white/[0.07] bg-[#0b1927]/90 px-4 py-[17px] shadow-[0_10px_35px_rgba(0,0,0,0.12)] transition duration-200 hover:border-white/[0.11] hover:bg-[#0d1c2b]">
      <div className="flex items-center gap-3.5">
        <div
          className={[
            "flex h-[43px] w-[43px] shrink-0 items-center justify-center rounded-full",
            iconBackgroundClassName,
          ].join(" ")}
        >
          <Icon
            size={22}
            strokeWidth={1.8}
            className={iconClassName}
          />
        </div>

        <div className="min-w-0">
          <p className="text-[11px] font-medium text-slate-400">
            {label}
          </p>

          <p
            className={[
              "mt-1 truncate text-[18px] font-semibold tracking-[-0.025em]",
              valueClassName,
            ].join(" ")}
            title={value}
          >
            {value}
          </p>

          <p className="mt-1 truncate text-[11px] text-slate-500">
            {description}
          </p>
        </div>
      </div>

      <div className="pointer-events-none absolute inset-x-0 top-0 h-px bg-gradient-to-r from-transparent via-white/[0.04] to-transparent" />
    </div>
  );
}
