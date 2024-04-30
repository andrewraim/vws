#ifndef REJECTION_CONTROL_H
#define REJECTION_CONTROL_H

namespace vws {

enum MaxRejectsAction { stop, warning, message, null };

//' Rejection Control
//'
//' Control object for rejection sampler.
//'
//' @param max_rejects Maximum number of rejections to tolerate before bailing out.
//' @param report Report progress each time this many candidates are proposed.
//' @param extra_outputs If \code{TRUE}, return a list with extended output
//' in addition to the accepted draws. Otherwise only return accepted draws.
//' @param action_incomplete What should happen if sampler halts with
//' \code{max_rejects} rejections: ne of \code{"stop"},  \code{"warning"}, or
//' \code{"message"}.
//'
//' @return
//' A control object to be passed to the \code{rejection} function.
//'
//' @name rejection_control
//' @export
class RejectionControl
{
public:
	RejectionControl(
		unsigned int max_rejects = 1000,
		unsigned int report_period = 100,
		MaxRejectsAction max_rejects_action = stop)
		: _max_rejects(max_rejects), _report_period(report_period),
		_max_rejects_action(max_rejects_action)
	{
	}

	unsigned int get_max_rejects() const { return _max_rejects; }
	unsigned int get_report_period() const { return _report_period; }
	MaxRejectsAction get_max_rejects_action() const { return _max_rejects_action; }

private:
	unsigned int _max_rejects;
	unsigned int _report_period;
	MaxRejectsAction _max_rejects_action;
};

}

#endif
