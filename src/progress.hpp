#ifndef PROGRESS_HPP
#define PROGRESS_HPP
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;

class ProgressBar
{
public:
    explicit ProgressBar(size_t total)
    {
        mTotal = total;
    }
    explicit ProgressBar(size_t total, bool display)
        : ProgressBar(total)
    {
        mIsDisplay = display;
    }

    void display()
    {
        if (mIsDisplay) Rcpp::Rcout << "| 20%| 40%| 60%| 80%| 100%|\n";
    }

    void tic()
    {
        mCurrent++;
        size_t percent = mCurrent * 100 / mTotal;
        if (mIsDisplay && percent > mPercentOld) {
            if (percent == 100) {
                Rcpp::Rcout << '|' << string(25, '*') << '|' << '\n';
            } else {
                size_t nbars = percent / 4;
                size_t vcur = percent % 4;
                string cur = " ";
                switch (vcur)
                {
                case 1:
                    cur = "•";
                    break;
                case 2:
                    cur = "-";
                    /* code */
                    break;
                case 3:
                    cur = "+";
                    /* code */
                    break;
                default:
                    break;
                }
                size_t nspaces = 25 - nbars - 1;
                string bar(nbars, '*'), space(nspaces, ' ');
                Rcpp::Rcout << '|' << bar << cur << space << '|' << string("\r");
            }
            mPercentOld = percent;
        }
    }

private:
    size_t mTotal;
    size_t mCurrent = 0;
    size_t mPercentOld = 0;
    bool mIsDisplay = true;
};

#endif  // PROGRESS_HPP