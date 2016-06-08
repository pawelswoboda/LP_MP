#ifndef LP_MP_ASCII_PLOT_VISITOR_HXX
#define LP_MP_ASCII_PLOT_VISITOR_HXX

#include "standard_visitor.hxx"
namespace ncurses {
#include <ncurses.h>
}
//#include <form.h>

namespace LP_MP {

// visitor deriving from StandardVisitor and printing convergence plots via curses and using ascii art
// also templatize for visitor class to make tightening visitor possible
template<class PROBLEM_DECOMPOSITION, class BASE_VISITOR = StandardVisitor<PROBLEM_DECOMPOSITION>>
class AsciiPlotVisitor : public BASE_VISITOR {
public:
   using BaseVisitor = BASE_VISITOR;
   // for collecting upper and lower bound per iteration
   struct Point {
      INDEX x;
      REAL y;
   };

   AsciiPlotVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd) 
      : BaseVisitor(cmd,pd)
   {}

   LPVisitorReturnType begin(const LP* lp)
   {
      auto ret = BaseVisitor::begin(lp);
      //const REAL lowerBound = StandardVisitor::GetLowerBound();
      //const REAL upperBound = lp->BestPrimalBound(); // note: LP does not compute initial primal bound right now. Implement this
      spdlog::get("logger")->info() << "note: LP does not compute initial primal bound right now. Implement this\n";
      ncurses::initscr(); 
      ncurses::start_color();
      ncurses::init_pair(1, COLOR_GREEN, COLOR_BLACK);
      ncurses::init_pair(2, COLOR_YELLOW, COLOR_BLACK);
      ncurses::clear();
      ncurses::refresh();

      view.ymin = BaseVisitor::GetLowerBound(); // do zrobienia: should be stored in LP
      view.ymax = lp->BestPrimalBound();  view.ymax = -30000;
      view.xmin = 0;
      view.xmax = 100;
      drawAxes(ncurses::stdscr,&view);
      ncurses::refresh();

      lowerBoundHistory_.push_back({BaseVisitor::GetIter(), BaseVisitor::GetLowerBound()});
      upperBoundHistory_.push_back({BaseVisitor::GetIter(), lp->BestPrimalBound()});

      return ret;
   }

   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(LP* lp)
   {
      auto ret_state = BaseVisitor::template visit<LP_STATE>(lp);
      const INDEX iter = BaseVisitor::GetIter();
      const REAL lowerBound = BaseVisitor::GetLowerBound();
      const REAL upperBound = lp->BestPrimalBound();
      lowerBoundHistory_.push_back({iter,lowerBound});
      upperBoundHistory_.push_back({iter,upperBound});

      // redraw graph with larger x-axis, if it does not fit into current plot
      if(iter >= view.xmax) {
         view.xmax *= 10;
         ncurses::erase();
         drawAxes(ncurses::stdscr, &view);
         drawGraph(ncurses::stdscr, &view, lowerBoundHistory_);
         drawGraph(ncurses::stdscr, &view, upperBoundHistory_);
         ncurses::refresh();
      }

      double d = 1;
      plotPoint(ncurses::stdscr, &view, REAL(BaseVisitor::GetIter()), lowerBound, '*', nullptr, nullptr);
      ncurses::refresh();
      return ret_state; 
   }

private:
   typedef double (*yfunction)(double x);
   struct _viewwin {
      double xmin, xmax;
      double ymin, ymax;
      double xscl, yscl;
   };
   typedef struct _viewwin viewwin;

   double estimateSlope(yfunction func, double x, double accuracy) const
   {
      double y1 = func(x - accuracy);
      double y2 = func(x + accuracy);
      return (y2 - y1) / (2 * accuracy);
   }
   double scale(double value, double omin, double omax, double nmin, double nmax)
   {
      /* Useful function to scale a value in one range to a different range.
         omin/omax - old range
         nmin/nmax - new range
         */
      double x = (value - omin) / (omax - omin);
      return x * (nmax - nmin) + nmin;
   }
   char slopeChar(double slope)
   {
      // Gets the character to display at a point in the graph with a given slope.
      double a = std::abs(slope);
      if (a < 0.5)        return '=';
      else if (a < 1.5)   return slope>0 ? '/' : '\\';
      else                return '|';
   }

   void plotPoint(ncurses::WINDOW *win, const viewwin *view, double x, double y, char ch, int *scrY, int *scrX)
   {
      /* Displays a point on the screen at a location determined by graph coordinates.
         win       - ncurses window for drawing (can be NULL to only set scrY and scrX w/o drawing)
         view      - view parameters structure
         x/y       - graph coordinates for point
         ch        - character to display
         scrY/scrX - screen coordinates where point was drawn are saved here if not NULL
         */
      int xm, ym; getmaxyx(win, ym, xm);
      int xp = scale(x, view->xmin, view->xmax, 0, xm);
      int yp = scale(y, view->ymin, view->ymax, ym, 0);
      spdlog::get("logger")->debug() << "print point (" << std::to_string(xp) << "," << std::to_string(yp) << ") with char " << ch ;

      if (scrX) *scrX = xp;
      if (scrY) *scrY = yp;

      assert(win);
      //if (win) 
      mvwaddch(win, yp, xp, ch);
      ncurses::refresh();
   }

   void getViewStep(ncurses::WINDOW *win, const viewwin *view, double *xstep, double *ystep)
   {
      // Gets the 'value' of one character on either or both axes.

      int xm, ym; getmaxyx(win, ym, xm);
      if (xstep) *xstep = (view->xmax - view->xmin) / (xm + 1);
      if (ystep) *ystep = (view->ymax - view->ymin) / (ym + 1);
   }

   void drawAxes(ncurses::WINDOW *win, const viewwin *view)
   {
      // This function is what draws the axes on the screen.

      int xm, ym; getmaxyx(win, ym, xm);
      double x0 = scale(0, view->xmin, view->xmax, 0, xm);
      double y0 = scale(0, view->ymin, view->ymax, ym, 0);

      double xstep, ystep; getViewStep(win, view, &xstep, &ystep);

      int i; for (i=0; i<=xm; i++) {
         double plotx = view->xmin + xstep * i;
         int tick = fabs(fmod(plotx, view->xscl)) < xstep;
         spdlog::get("logger")->info() << "print vertical point (" << std::to_string(y0) << "," << std::to_string(i) << ")";
         mvwaddch(win, y0, i, tick ? '+':'-');
      }
      for (i=0; i<=ym; i++) {
         double ploty = view->ymin + ystep * i;
         int tick = fabs(fmod(ploty, view->yscl)) < ystep;
         spdlog::get("logger")->info() << "print horizontal point (" << std::to_string(i) << "," << std::to_string(x0) << ")";
         mvwaddch(win, i, x0, tick ? '+':'|');
      }

      mvwaddch(win, y0, x0, '+');
   }

   void drawGraph(ncurses::WINDOW *win, const viewwin *view, const std::vector<Point>& graph)
   {
      /* Draws a graph on the screen without axes.
         win              - ncurses window for drawing
         view             - view parameters structure
         yfunc            - function to graph (function pointer)
         enableSlopeChars - whether or not to call slopeChar to determine characters
         */
      int xm, ym; getmaxyx(win, ym, xm);
      for(INDEX i=0; i<graph.size(); ++i) {
         double y = graph[i].y;
         //double d = estimateSlope(yfunc, x, step/2);
         plotPoint(win, view, graph[i].x, y, '*', NULL, NULL);
      }
   }






   double xMin, xMax;
   double yMin, yMax;
   viewwin view; //do zrobienia: make better
   // history of upper and lower bound, needed, when we redraw the window and for choosing slopes
   std::vector<Point> lowerBoundHistory_;
   std::vector<Point> upperBoundHistory_;
};

} // end namespace LP_MP

#endif // LP_MP_ASCII_PLOT_VISITOR_HXX
