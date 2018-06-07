#ifndef __GEODA_CENTER_GWR_DLG_H__
#define __GEODA_CENTER_GWR_DLG_H__

#include <vector>
#include <wx/dialog.h>
#include <wx/gauge.h>
#include <Eigen/Dense>

#include "../FramesManagerObserver.h"
#include "../DataViewer/TableStateObserver.h"
#include "../ShapeOperations/WeightsManStateObserver.h"
#include "RegressionReportDlg.h"
#include "SaveToTableDlg.h"


class FramesManager;
class TableState;
class TableInterface;
class Project;
struct SaveToTableEntry;
using namespace Eigen;




class GWRDlg:public wxDialog,public FramesManagerObserver,public TableStateObserver
{
	DECLARE_EVENT_TABLE()
public:
	GWRDlg(Project* project,
			wxWindow* parent=NULL,
			wxString title=_("GWR"),
			wxWindowID id=-1,
			const wxString& caption = _("GWR"),
			const wxPoint& pos=wxDefaultPosition,
			const wxSize& size=wxDefaultSize,
			long style = wxCAPTION|wxDEFAULT_DIALOG_STYLE);
    virtual ~GWRDlg();

	bool Create( wxWindow* parent, wxWindowID id = -1,
				const wxString& caption = _("GWR"),
				const wxPoint& pos = wxDefaultPosition,
				const wxSize& size = wxDefaultSize,
				long style = wxCAPTION|wxDEFAULT_DIALOG_STYLE );
	void CreateControls();
	//button event
	void OnApplyClick(wxCommandEvent& event);
	void OnCloseClick( wxCommandEvent& event );
	void OnClose(wxCloseEvent& event);
	void OnBtnLeft1( wxCommandEvent& event );
	void OnBtnRight1( wxCommandEvent& event );
	void OnBtnLeft2( wxCommandEvent& event );
	void OnBtnRight2( wxCommandEvent& event );
	void OnBtnLeft3( wxCommandEvent& event );
	void OnBtnRight3( wxCommandEvent& event );
	void OnBtnLeft4( wxCommandEvent& event );
	void OnBtnRight4( wxCommandEvent& event );
	void OnBtnLeft5( wxCommandEvent& event );
	void OnBtnRight5( wxCommandEvent& event );
	void OnRadioGaussianSel( wxCommandEvent& event );
	void OnRadioBisquareSel( wxCommandEvent& event );
	void OnRadioAICSel( wxCommandEvent& event );
	void OnRadioCVSel( wxCommandEvent& event );
	void OnbwtextChanged( wxCommandEvent& event );
	void ConstructDepMatrix();
	void ConstructCoorMatrix();
	void ConstructIndepMatrix();
	//bool OutCsv(MatrixXd outmat);
	//bool OutCsv(VectorXd outmat);

	VectorXd eu_dist_vec(MatrixXd in_locs,MatrixXd out_loc);
	VectorXd bisq_wt_vec(VectorXd distv, double bw);
	VectorXd gauss_wt_vec(VectorXd distv, double bw);
	MatrixXd Parameters_matrix(MatrixXd in_locs,double bw,int indepvarnum);
	VectorXd gw_reg(MatrixXd x,VectorXd y,MatrixXd w,bool hatmatrix,int focus);
	VectorXd ehat(VectorXd y, MatrixXd X, MatrixXd beta);
	double rss(VectorXd y, MatrixXd X, MatrixXd beta);
	VectorXd trhat2(MatrixXd S);
	std::vector<double> gwr_diag(VectorXd y,MatrixXd x, MatrixXd beta, MatrixXd S);
	bool CheckParameter();
	double Goldensection_search(double xlower,double xupper);
	VectorXd GetDistanceVec(MatrixXd in_locs);
	double ForGss();
	double GetAIC(double bw);
	double GetCV(double bw);
	void PrintGWRresult();//计算结果
	void DisplayRegression(wxString dump);//新窗口显示结果
	void SaveDataToTable();//保存beta进表
	void OnReportClose(wxWindowDestroyEvent& event);
	double Percentile(double x, const std::vector<double>& v);
	void BreaksPoint(const std::vector<double>& var);

	void InitVariableList();
	void EnablingItems();
	void UpdateMessageBox(wxString msg);

	/** Implementation of FramesManagerObserver interface */
	virtual void update(FramesManager* o) ;
	/** Implementation of TableStateObserver interface */
	virtual void update(TableState* o);
	virtual bool AllowTimelineChanges() { return true; }
	virtual bool AllowGroupModify(const wxString& grp_nm) { return true; }
	virtual bool AllowObservationAddDelete() { return false; }

	wxListBox* m_varlist;
	wxTextCtrl* m_dependent;
	wxTextCtrl* m_ID;
	wxTextCtrl* m_X;
	wxTextCtrl* m_Y;
	wxListBox* m_independentlist;
	wxRadioButton* m_gaussian;
	wxRadioButton* m_bisquare;
	wxRadioButton* m_aic;
	wxRadioButton* m_cv;
	wxTextCtrl* m_bandwidth;
	RegressionReportDlg *regReportDlg;
	wxGauge* m_gauge;
	wxStaticText* m_gauge_text;
	long		m_obs;
	int			WeightModel;
	int			BandWidthModel;
	double		bw;
	bool		m_apply;
	std::map<wxString, wxString> name_to_nm;
	std::map<wxString, int> name_to_tm_id;
	std::vector<double> diagnosticVal;
	wxString logReport;
	std::vector<double> breaks;

	VectorXd dep_vec;
	MatrixXd indep_mat;
	MatrixXd coor_mat;
	MatrixXd S;
	MatrixXd yhat;
	MatrixXd beta_se;
	MatrixXd outbetamat;

	Project* project;
	TableInterface* table_int;
	FramesManager* frames_manager;

};

#endif
