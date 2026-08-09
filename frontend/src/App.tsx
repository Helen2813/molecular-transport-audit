import {
  useState,
} from "react";

import NewRunModal from "./components/dialogs/NewRunModal";
import Sidebar from "./components/layout/Sidebar";
import TopBar from "./components/layout/TopBar";
import OverviewPage from "./pages/OverviewPage";


export default function App() {
  const [
    newRunOpen,
    setNewRunOpen,
  ] = useState(
    false,
  );

  return (
    <div className="min-h-screen bg-transparent text-white">
      <Sidebar />

      <TopBar
        onNewRun={() =>
          setNewRunOpen(
            true,
          )
        }
      />

      <OverviewPage />

      <NewRunModal
        open={
          newRunOpen
        }
        onClose={() =>
          setNewRunOpen(
            false,
          )
        }
      />
    </div>
  );
}
