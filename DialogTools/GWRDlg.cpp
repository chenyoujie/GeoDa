#include <time.h>
#include <boost/foreach.hpp>
#include <wx/wx.h>
#include <wx/grid.h>
#include <wx/msgdlg.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/xrc/xmlres.h> // XRC XML resouces
#include <wx/filedlg.h>
#include <wx/textdlg.h>
#include <iostream>
#include <fstream>
#include <math.h>


#include "../FramesManager.h"
#include "../GenUtils.h"
#include "../GeoDa.h"
#include "../Project.h"
#include "../TemplateCanvas.h"
#include "GWRDlg.h"
#include "SaveToTableDlg.h"
#include "../DataViewer/TableInterface.h"
#include "../DataViewer/TableState.h"
#include "../DbfFile.h"
#include <wx/msgdlg.h>
#include "RegressionReportDlg.h"


BEGIN_EVENT_TABLE( GWRDlg, wxDialog )
	EVT_BUTTON(XRCID("ID_BTN_CANCEL"),GWRDlg::OnCloseClick)
	EVT_BUTTON(XRCID("IDG_BTN_LEFT1"),GWRDlg::OnBtnLeft1)
	EVT_BUTTON(XRCID("IDG_BTN_RIGHT1"),GWRDlg::OnBtnRight1)
	EVT_BUTTON(XRCID("IDG_BTN_LEFT2"),GWRDlg::OnBtnLeft2)
	EVT_BUTTON(XRCID("IDG_BTN_RIGHT2"),GWRDlg::OnBtnRight2)
	EVT_BUTTON(XRCID("IDG_BTN_LEFT3"),GWRDlg::OnBtnLeft3)
	EVT_BUTTON(XRCID("IDG_BTN_RIGHT3"),GWRDlg::OnBtnRight3)
	EVT_BUTTON(XRCID("IDG_BTN_LEFT4"),GWRDlg::OnBtnLeft4)
	EVT_BUTTON(XRCID("IDG_BTN_RIGHT4"),GWRDlg::OnBtnRight4)
	EVT_BUTTON(XRCID("IDG_BTN_LEFT5"),GWRDlg::OnBtnLeft5)
	EVT_BUTTON(XRCID("IDG_BTN_RIGHT5"),GWRDlg::OnBtnRight5)
	EVT_RADIOBUTTON(XRCID("IDG_RADBTN_GAUSSIAN"),GWRDlg::OnRadioGaussianSel)
	EVT_RADIOBUTTON(XRCID("IDG_RADBTN_BISQ"),GWRDlg::OnRadioBisquareSel)
	EVT_RADIOBUTTON(XRCID("IDG_RADBTN_AIC"),GWRDlg::OnRadioAICSel)
	EVT_RADIOBUTTON(XRCID("IDG_RADBTN_CV"),GWRDlg::OnRadioCVSel)
	EVT_BUTTON(XRCID("IDG_BTN_APPLY"),GWRDlg::OnApplyClick)
	EVT_TEXT(XRCID("IDG_TEXTCTRL_BW"),GWRDlg::OnbwtextChanged)
END_EVENT_TABLE()

GWRDlg::GWRDlg(Project* project_s,wxWindow* parent ,
				wxString title,
				wxWindowID id,const wxString& caption,
				const wxPoint& pos,const wxSize& size,
				long style)
: project(project_s), frames_manager(project_s->GetFramesManager()),
table_int(project_s->GetTableInt()),
regReportDlg(0)

{
	wxLogMessage("open GWRDlg.");
	//project = project_s;
	Create(parent,id,caption,pos,size,style);
	InitVariableList();
}
GWRDlg::~GWRDlg()
{
	wxLogMessage("GWRDlg::~GWRDlg()");
}

bool GWRDlg::Create(wxWindow* parent, wxWindowID id,
					const wxString& caption, const wxPoint& pos,
					const wxSize& size, long style )
{
	m_varlist=NULL;
	m_dependent = NULL;
	m_ID = NULL;
	m_X = NULL;
	m_Y = NULL;;
	m_independentlist = NULL;
	m_gaussian = NULL;
	m_bisquare = NULL;
	m_aic = NULL;
	m_cv = NULL;
	m_bandwidth = NULL;
	//m_gauge = NULL;
	m_gauge_text = NULL;
	WeightModel=0;
	BandWidthModel=0;
	SetParent(parent);
	CreateControls();
	Centre();

	return true;
}

void GWRDlg::CreateControls()
{    
	wxXmlResource::Get()->LoadDialog(this, GetParent(), "IDG_GWRMODEL");
	m_varlist = XRCCTRL(*this, "IDG_LIST_VAR", wxListBox);
	m_dependent = XRCCTRL(*this, "IDG_DEPEND_VAR" ,wxTextCtrl);
	m_ID = XRCCTRL(*this, "IDG_ID" ,wxTextCtrl);
	m_X = XRCCTRL(*this, "IDG_COR_X" ,wxTextCtrl);
	m_Y = XRCCTRL(*this, "IDG_COR_Y" ,wxTextCtrl);
	m_independentlist = XRCCTRL(*this, "IDG_INDEP_VAR", wxListBox);
	m_gaussian = XRCCTRL(*this, "IDG_RADBTN_GAUSSIAN", wxRadioButton);
	m_bisquare = XRCCTRL(*this, "IDG_RADBTN_BISQ", wxRadioButton);
	m_aic = XRCCTRL(*this, "IDG_RADBTN_AIC", wxRadioButton);
	m_cv = XRCCTRL(*this, "IDG_RADBTN_CV", wxRadioButton);
	m_bandwidth = XRCCTRL(*this, "IDG_TEXTCTRL_BW" ,wxTextCtrl);
	m_gauge = XRCCTRL(*this, "IDG_GAUGE", wxGauge);
	m_gauge->SetRange(200);
	m_gauge->SetValue(50);
	m_gauge->Hide();
	m_gauge_text = XRCCTRL(*this, "IDG_GAUGE_TEXT", wxStaticText);
}

void GWRDlg::InitVariableList()
{
	wxLogMessage("GWRDlg::InitVariableList()");
	UpdateMessageBox(wxEmptyString);

	m_obs = project->GetNumRecords();
	if (m_obs <= 0) {
		wxMessageBox("Error: no records found in DBF file");
		return;
	}

	std::vector<int> col_id_map;
	table_int->FillNumericColIdMap(col_id_map);
	name_to_nm.clear(); // map to table_int col id
	name_to_tm_id.clear(); // map to corresponding time id
	for (int i=0, iend=col_id_map.size(); i<iend; i++) {
		int id = col_id_map[i];
		wxString name = table_int->GetColName(id);
		if (table_int->IsColTimeVariant(id)) {
			for (int t=0; t<table_int->GetColTimeSteps(id); t++) {
				wxString nm = name;
				nm << " (" << table_int->GetTimeString(t) << ")";
				name_to_nm[nm] = name;
				name_to_tm_id[nm] = t;
				m_varlist->Append(nm);
			}
		} 
		else {
			name_to_nm[name] = name;
			name_to_tm_id[name] = 0;
			m_varlist->Append(name);
		}
	}

	//y = NULL;
	//x = NULL;
	m_varlist->SetSelection(0);
	m_varlist->SetFirstItem(m_varlist->GetSelection());
	m_apply=false;
	EnablingItems();
}

void GWRDlg::UpdateMessageBox(wxString msg)
{
	m_gauge_text->SetLabelText(msg);
	m_gauge_text->Update();
}

void GWRDlg::EnablingItems()
{
	FindWindow(XRCID("IDG_BTN_APPLY"))->Enable(m_apply);	
	if (m_ID->GetLineText(0)==wxEmptyString)
	{
		FindWindow(XRCID("IDG_BTN_LEFT1"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT1"))->Enable(false);
	}
	else
	{
		FindWindow(XRCID("IDG_BTN_LEFT1"))->Enable(false);
		FindWindow(XRCID("IDG_BTN_RIGHT1"))->Enable(true);
	}
	if (m_X->GetLineText(0)==wxEmptyString)
	{
		FindWindow(XRCID("IDG_BTN_LEFT2"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT2"))->Enable(false);
	}
	else
	{
		FindWindow(XRCID("IDG_BTN_LEFT2"))->Enable(false);
		FindWindow(XRCID("IDG_BTN_RIGHT2"))->Enable(true);
	}
	if (m_Y->GetLineText(0)==wxEmptyString)
	{
		FindWindow(XRCID("IDG_BTN_LEFT3"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT3"))->Enable(false);
	}
	else
	{
		FindWindow(XRCID("IDG_BTN_LEFT3"))->Enable(false);
		FindWindow(XRCID("IDG_BTN_RIGHT3"))->Enable(true);
	}
	if (m_dependent->GetLineText(0)==wxEmptyString)
	{
		FindWindow(XRCID("IDG_BTN_LEFT4"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT4"))->Enable(false);
	}
	else
	{
		FindWindow(XRCID("IDG_BTN_LEFT4"))->Enable(false);
		FindWindow(XRCID("IDG_BTN_RIGHT4"))->Enable(true);
	}
	if (m_independentlist->GetCount()==0)
	{
		FindWindow(XRCID("IDG_BTN_LEFT5"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT5"))->Enable(false);
	}
	else
	{
		FindWindow(XRCID("IDG_BTN_LEFT5"))->Enable(true);
		FindWindow(XRCID("IDG_BTN_RIGHT5"))->Enable(true);
	}
}
void GWRDlg::update(FramesManager* o)
{
	wxLogMessage("update");
}
void GWRDlg::update(TableState* o)
{
	wxLogMessage("update");
}

void GWRDlg::OnCloseClick( wxCommandEvent& event )
{
	wxLogMessage("Click GWRDlg::OnCloseClick");

	event.Skip();
	EndDialog(wxID_CLOSE);
	Destroy();
}

void GWRDlg::OnClose(wxCloseEvent& event)
{
	Destroy();
}

void GWRDlg::OnBtnLeft1( wxCommandEvent& event )
{
	if(m_varlist->GetCount() > 0)
	{
		if (m_varlist ->GetSelection() >=0)
		{
			wxString temp=m_ID ->GetValue();
			m_ID->SetValue(m_varlist->GetString(m_varlist->GetSelection()));
			m_varlist->Delete(m_varlist->GetSelection());
		}
	}
	EnablingItems();
}
void GWRDlg::OnBtnRight1( wxCommandEvent& event )
{
	if (m_ID->GetLineText(0)!=wxEmptyString)
	{
		m_varlist->Append(m_ID->GetLineText(0));
		m_ID->Clear();
	}
	EnablingItems();
}

void GWRDlg::OnBtnLeft2( wxCommandEvent& event )
{
	if(m_varlist->GetCount() > 0)
	{
		if (m_varlist ->GetSelection() >=0)
		{
			wxString temp=m_X ->GetValue();
			m_X->SetValue(m_varlist->GetString(m_varlist->GetSelection()));
			m_varlist->Delete(m_varlist->GetSelection());
		}
	}
	EnablingItems();
}
void GWRDlg::OnBtnRight2( wxCommandEvent& event )
{
	if (m_X->GetLineText(0)!=wxEmptyString)
	{
		m_varlist->Append(m_X->GetLineText(0));
		m_X->Clear();
	}
	EnablingItems();
}

void GWRDlg::OnBtnLeft3( wxCommandEvent& event )
{
	if(m_varlist->GetCount() > 0)
	{
		if (m_varlist ->GetSelection() >=0)
		{
			wxString temp=m_Y ->GetValue();
			m_Y->SetValue(m_varlist->GetString(m_varlist->GetSelection()));
			m_varlist->Delete(m_varlist->GetSelection());
		}
	}
	EnablingItems();
}
void GWRDlg::OnBtnRight3( wxCommandEvent& event )
{
	if (m_Y->GetLineText(0)!=wxEmptyString)
	{
		m_varlist->Append(m_Y->GetLineText(0));
		m_Y->Clear();
	}
	EnablingItems();
}

void GWRDlg::OnBtnLeft4( wxCommandEvent& event )
{
	if(m_varlist->GetCount() > 0)
	{
		if (m_varlist ->GetSelection() >=0)
		{
			wxString temp=m_dependent ->GetValue();
			m_dependent->SetValue(m_varlist->GetString(m_varlist->GetSelection()));
			m_varlist->Delete(m_varlist->GetSelection());
		}
	}
	EnablingItems();
}
void GWRDlg::OnBtnRight4( wxCommandEvent& event )
{
	if (m_dependent->GetLineText(0)!=wxEmptyString)
	{
		m_varlist->Append(m_dependent->GetLineText(0));
		m_dependent->Clear();
	}
	EnablingItems();
}

void GWRDlg::OnBtnLeft5( wxCommandEvent& event )
{
	if(m_varlist->GetCount() > 0)
	{
		if (m_varlist ->GetSelection() >=0)
		{
			m_independentlist->Append(m_varlist->GetString(m_varlist->GetSelection()));
			m_varlist->Delete(m_varlist->GetSelection());
		}
	}
	EnablingItems();
}
void GWRDlg::OnBtnRight5( wxCommandEvent& event )
{
	if (m_independentlist->GetCount()>0)
	{
		if(m_independentlist->GetSelection()>=0)
		{
			int cur_sel=m_independentlist->GetSelection();
			m_varlist->Append(m_independentlist->GetString(cur_sel));
			m_independentlist->Delete(cur_sel);
		}
	}
	EnablingItems();
}

void GWRDlg::OnRadioGaussianSel( wxCommandEvent& event )
{
	WeightModel=1;
}
void GWRDlg::OnRadioBisquareSel( wxCommandEvent& event )
{
	WeightModel=2;
}
void GWRDlg::OnRadioAICSel( wxCommandEvent& event )
{
	BandWidthModel=1;
	FindWindow(XRCID("IDG_TEXTCTRL_BW"))->Enable(false);
	m_apply=CheckParameter();
	EnablingItems();
}
void GWRDlg::OnRadioCVSel( wxCommandEvent& event )
{
	BandWidthModel=2;
	FindWindow(XRCID("IDG_TEXTCTRL_BW"))->Enable(false);
	m_apply=CheckParameter();
	EnablingItems();
}
void GWRDlg::OnbwtextChanged( wxCommandEvent& event )
{
	m_apply=CheckParameter();
	EnablingItems();
}

void GWRDlg::OnApplyClick(wxCommandEvent& event)
{
	ConstructIndepMatrix();
	ConstructCoorMatrix();
	ConstructDepMatrix();
	if (CheckParameter())
	{
		m_gauge->Show();
		m_gauge->Update();
		UpdateMessageBox("calculating...");
		if (BandWidthModel==0)
		{
			MatrixXd beta_matrix=Parameters_matrix(coor_mat,bw,indep_mat.cols());
			vector<double> result1=gwr_diag(dep_vec,indep_mat,beta_matrix,S);
		}
		else 
			bw=ForGss();
		m_gauge->SetValue(150);
		m_gauge->Update();
		SaveDataToTable();
		PrintGWRresult();
		DisplayRegression(logReport);
		m_gauge->SetValue(200);
		UpdateMessageBox("done");
	}
}

void GWRDlg::ConstructIndepMatrix()
{
	const unsigned int indepnum = m_independentlist ->GetCount();//生成自变量矩阵
	std::vector<double> vec(m_obs);
	indep_mat.resize(m_obs,indepnum+1);
	for (int i=0;i<indepnum;i++) 
	{
		wxString nm=name_to_nm[m_independentlist->GetString(i)];
		int col=table_int->FindColId(nm);
		if (col==wxNOT_FOUND)
		{
			wxString err_msg = wxString::Format(_("Variable %s is no longer in the Table.  Please close and reopen the Regression Dialog to synchronize with Table data."), nm);
			wxMessageDialog dlg(NULL, err_msg, "Error", wxOK | wxICON_ERROR);
			dlg.ShowModal();
			return;
		}
		int tm=name_to_tm_id[m_independentlist->GetString(i)];
		table_int->GetColData(col,tm,vec);
		if (i==0)
		{
			for(int l=0;l<m_obs;l++)
				indep_mat(l,i)=1;
		}
		for (int j=0;j<m_obs;j++)
		{
			indep_mat(j,i+1)=vec[j];
		}
	}
}
void GWRDlg::ConstructDepMatrix()
{
	std::vector<double> vec(m_obs);
	dep_vec.resize(m_obs);
	wxString nm=name_to_nm[m_dependent->GetLineText(0)];
	int col=table_int->FindColId(nm);
	if (col==wxNOT_FOUND)
	{
		wxString err_msg = wxString::Format(_("Variable %s is no longer in the Table.  Please close and reopen the Regression Dialog to synchronize with Table data."), nm);
		wxMessageDialog dlg(NULL, err_msg, "Error", wxOK | wxICON_ERROR);
		dlg.ShowModal();
		return;
	}
	int tm=name_to_tm_id[m_dependent->GetLineText(0)];
	table_int->GetColData(col,tm,vec);
	for (int i=0;i<m_obs;i++)
	{
		dep_vec(i)=vec[i];
	}
	/*if (OutCsv(dep_vec))		wxMessageBox("success!");*/
}
void GWRDlg::ConstructCoorMatrix()
{
	std::vector<double> vecX(m_obs);//生成坐标矩阵
	std::vector<double> vecY(m_obs);
	coor_mat.resize(m_obs,2);
	wxString nmX=name_to_nm[m_X->GetLineText(0)];
	int colX=table_int->FindColId(nmX);
	wxString nmY=name_to_nm[m_Y->GetLineText(0)];
	int colY=table_int->FindColId(nmY);
	if (colX==wxNOT_FOUND)
	{
		wxString err_msg = wxString::Format(_("Variable %s is no longer in the Table.  Please close and reopen the Regression Dialog to synchronize with Table data."), nmX);
		wxMessageDialog dlg(NULL, err_msg, "Error", wxOK | wxICON_ERROR);
		dlg.ShowModal();
		return;
	}
	if (colY==wxNOT_FOUND)
	{		
		wxString err_msg = wxString::Format(_("Variable %s is no longer in the Table.  Please close and reopen the Regression Dialog to synchronize with Table data."), nmY);
		wxMessageDialog dlg(NULL, err_msg, "Error", wxOK | wxICON_ERROR);
		dlg.ShowModal();
		return;
	}
	int tmX=name_to_tm_id[m_X->GetLineText(0)];
	table_int->GetColData(colX,tmX,vecX);
	int tmY=name_to_tm_id[m_Y->GetLineText(0)];
	table_int->GetColData(colY,tmY,vecY);
	for (int i=0;i<m_obs;i++)
	{
		coor_mat(i,0)=vecX[i];
		coor_mat(i,1)=vecY[i];
	}	
	
}
bool GWRDlg::CheckParameter()
{
	bool result=true;
	if (m_ID->GetLineText(0)==wxEmptyString)
	{wxMessageBox("ID is necessary");result=false;}
	if (m_X->GetLineText(0)==wxEmptyString)
	{wxMessageBox("X is necessary");result=false;}
	if (m_Y->GetLineText(0)==wxEmptyString)
	{wxMessageBox("Y is necessary");result=false;}
	if (m_dependent->GetLineText(0)==wxEmptyString)
	{wxMessageBox("dependent variable is necessary");result=false;}
	if (m_independentlist->GetCount()==0)
	{wxMessageBox("independent variable is necessary");result=false;}
	if (WeightModel==0)
	{wxMessageBox("please choose a weight model!");result=false;}
	if (m_bandwidth->GetLineText(0)==wxEmptyString&&BandWidthModel==0)
	{wxMessageBox("bandwidth is necessary");result=false;}
	wxString str_bw=m_bandwidth->GetLineText(0);
	if (str_bw!=wxEmptyString)
		str_bw.ToDouble(&bw);
	return result;
}

VectorXd GWRDlg::eu_dist_vec(MatrixXd in_locs,MatrixXd out_loc)
{
	int n_in = in_locs.rows();
	VectorXd eu_dist(n_in);
	double inner;
	for(int i=0;i<n_in;i++){
		inner=(in_locs.row(i).array()-out_loc.array()).pow(2).sum();
		eu_dist(i) = sqrt(inner);
	}
	return eu_dist;
}

VectorXd GWRDlg::bisq_wt_vec(VectorXd distv, double bw)
{
	int n = distv.size();
	VectorXd wtv;
	wtv.setZero(n);
	for(int i=0; i<n; i++){
		if(distv(i) <= bw)
			wtv(i) = pow(1 - pow(distv(i)/bw,2),2);
	}
	return wtv;
}

VectorXd GWRDlg::gauss_wt_vec(VectorXd distv, double bw)
{
	int n = distv.size();
	VectorXd wtv;
	wtv.setZero(n);
	for(int i=0; i<n; i++){
		wtv(i) = exp(pow(distv(i),2)/((-2)*pow(bw,2)));
	}
	return wtv;
}

VectorXd GWRDlg::gw_reg(MatrixXd x,VectorXd y,MatrixXd w,bool hatmatrix,int focus)
{
	MatrixXd xtw=x.transpose()*w;
	VectorXd beta=(xtw*x).colPivHouseholderQr().solve(xtw*y);

	if(hatmatrix)
	{
		if(focus==0)
		{
			MatrixXd xtwx_inv=(xtw*x).inverse();
			MatrixXd ci=xtwx_inv*xtw;
			MatrixXd s_ri =x.row(focus)*ci;
			VectorXd beta_sev = (ci*(ci.transpose())).diagonal();
			beta_se.resize(y.size(),beta_sev.size());
			beta_se.row(focus)=beta_sev.transpose();
			S.resize(x.rows(),ci.cols());
			S.row(focus)=s_ri;
		}
		else
		{
			MatrixXd xtwx_inv=(xtw*x).inverse();
			MatrixXd ci=xtwx_inv*xtw;
			MatrixXd s_ri =x.row(focus)*ci;
			VectorXd beta_sev = (ci*(ci.transpose())).diagonal();
			beta_se.row(focus)=beta_sev.transpose();
			S.row(focus)=s_ri;
		}
	}
	return beta;
}


MatrixXd GWRDlg::Parameters_matrix(MatrixXd in_locs,double bw,int indepvarnum)
{
	int sampleNum=in_locs.rows();
	MatrixXd wt_matrix;
	VectorXd eu_dist;
	MatrixXd out_loc;
	VectorXd weight_mat;
	VectorXd outbetavec;
	outbetamat.resize(sampleNum,indepvarnum);
	if (BandWidthModel==2)
	{
		for(int i=0;i<sampleNum;i++)
		{
			out_loc=in_locs.row(i);
			eu_dist=eu_dist_vec(in_locs,out_loc);
			switch(WeightModel)
			{
			case 1:
				weight_mat=gauss_wt_vec(eu_dist,bw);
				break;
			case 2:
				weight_mat=bisq_wt_vec(eu_dist,bw);
				break;
			}
			weight_mat(i)=0;
			wt_matrix=weight_mat.asDiagonal();
			outbetavec=gw_reg(indep_mat,dep_vec,wt_matrix,true,i);
			outbetamat.row(i)=outbetavec.transpose();
		}
		yhat=S*dep_vec;
	}
	else
	{
		for(int i=0;i<sampleNum;i++)
		{
			out_loc=in_locs.row(i);
			eu_dist=eu_dist_vec(in_locs,out_loc);
			switch(WeightModel)
			{
			case 1:
				weight_mat=gauss_wt_vec(eu_dist,bw);
				break;
			case 2:
				weight_mat=bisq_wt_vec(eu_dist,bw);
				break;
			}
			wt_matrix=weight_mat.asDiagonal();
			outbetavec=gw_reg(indep_mat,dep_vec,wt_matrix,true,i);
			outbetamat.row(i)=outbetavec.transpose();
		}
		yhat=S*dep_vec;
	}
	return outbetamat;
}

VectorXd GWRDlg::ehat(VectorXd y, MatrixXd X, MatrixXd beta) 
{
	VectorXd fitted = (beta.cwiseProduct(X)).rowwise().sum();
	return y - fitted;
}

double GWRDlg::rss(VectorXd y, MatrixXd X, MatrixXd beta) 
{
	VectorXd r = ehat(y, X, beta);
	return (r.cwiseProduct(r)).sum();
}

VectorXd GWRDlg::trhat2(MatrixXd S) 
{
	int n_obs = S.rows();
	double htr = 0.0;
	double htr2 = 0.0;
	VectorXd result(2);
	for (int i = 0; i < n_obs; i++) {
		htr  += S(i,i);
		htr2 += (S.row(i).cwiseProduct(S.row(i))).sum();
	}
	result(0) = htr;
	result(1) = htr2;
	return result;
}

vector<double> GWRDlg::gwr_diag(VectorXd y,MatrixXd x, MatrixXd beta, MatrixXd S)
{
	double ss = rss(y,x,beta);
	VectorXd s_hat = trhat2(S);
	int n = S.rows();
	/*vector<double> result;*/
	double AIC = n*log(ss/n)+n*log(2*M_PI)+n+s_hat(0); //AIC
	double AICc = n*log(ss/n)+n*log(2*M_PI)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
	double edf = n-2*s_hat(0) + s_hat(1); //edf
	double enp = 2*s_hat(0) - s_hat(1); // enp
	VectorXd meany(n);
	meany.fill(y.sum()/n);
	double yss = ((y.array()-meany.array()).pow(2)).sum(); //yss.g
	double r2 = 1 - ss/yss;
	double r2_adj = 1-(1-r2)*(n-1)/(edf-1);
	diagnosticVal.clear();
	diagnosticVal.push_back(enp) ;
	diagnosticVal.push_back(s_hat(0)) ;
	diagnosticVal.push_back(s_hat(1)) ;

	diagnosticVal.push_back(r2) ;
	diagnosticVal.push_back(r2_adj) ;
	diagnosticVal.push_back(AIC) ;
	diagnosticVal.push_back(AICc) ;
	diagnosticVal.push_back(edf) ;
	diagnosticVal.push_back(ss) ;

	return diagnosticVal;
}

VectorXd GWRDlg::GetDistanceVec(MatrixXd in_locs)
{
	int sampleNum=in_locs.rows();
	VectorXd eu_dist;
	MatrixXd out_loc;
	out_loc=in_locs.row(0);
	eu_dist=eu_dist_vec(in_locs,out_loc);
	return eu_dist;
}
double GWRDlg::GetCV(double bw)
{
	int indepvarnum=indep_mat.cols();
	MatrixXd beta_matrix=Parameters_matrix(coor_mat,bw,indepvarnum);
	//vector<double> result = gwr_diag(dep_vec,indep_mat,beta_matrix,S);
	double CV;
	CV=((dep_vec.array()-yhat.array()).pow(2)).sum();
	return CV;
}
double GWRDlg::GetAIC(double bw)
{
	int indepvarnum=indep_mat.cols();
	double AICc;
	MatrixXd beta_matrix=Parameters_matrix(coor_mat,bw,indepvarnum);
	//vector<double> result = gwr_diag(dep_vec,indep_mat,beta_matrix,S);
	//AICc=result[6];
	double ss = rss(dep_vec,indep_mat,beta_matrix);
	VectorXd s_hat = trhat2(S);
	int n = S.rows();
	AICc = n*log(ss/n)+n*log(2*M_PI)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc

	return AICc;
}

double GWRDlg::Goldensection_search(double xlower,double xupper)
{
	double eps=1e-4;
	double a=5;
	double R=(sqrt(a)-1)/2;
	double xopt;
	double x1,x2;
	double f1,f2;
	double d=R*(xupper-xlower);
	x1=xlower+d;
	x2=xupper-d;
	if(BandWidthModel==1)
	{
		f1=GetAIC(x1);
		f2=GetAIC(x2);
		double d1=f2-f1;
		if(f1<f2)
			xopt=x1;
		else
			xopt=x2;
		while(abs(d)>eps && abs(d1)>eps)  //AIC的函数，是单调的吗？
		{
			d=R*d;
			if(f1<f2)
			{
				xlower=x2;
				x2=x1;
				x1=xlower+d;
				f2=f1;
				f1=GetAIC(x1);
			}
			else
			{
				xupper=x1;
				x1=x2;
				x2=xupper-d;
				f1=f2;
				f2=GetAIC(x2);
			}
		}
		if(f1<f2)
			xopt=x1;
		else
			xopt=x2;
	}
	if(BandWidthModel==2)
	{
		f1=GetCV(x1);
		f2=GetCV(x2);
		double d1=f2-f1;
		if(f1<f2)
			xopt=x1;
		else
			xopt=x2;
		while(abs(d)>eps && abs(d1)>eps)
		{
			d=R*d;
			if(f1<f2)
			{
				xlower=x2;
				x2=x1;
				x1=xlower+d;
				f2=f1;
				f1=GetCV(x1);
			}
			else
			{
				xupper=x1;
				x1=x2;
				x2=xupper-d;
				f1=f2;
				f2=GetCV(x2);
			}
		}
		if(f1<f2)
			xopt=x1;
		else
			xopt=x2;
	}
	return xopt;
}

double GWRDlg::ForGss()
{
	VectorXd eu_dist=GetDistanceVec(coor_mat);
	double upper = eu_dist.maxCoeff();
	double lower = upper/5000;
	double xopt = Goldensection_search(lower,upper);
	vector<double> result = gwr_diag(dep_vec,indep_mat,outbetamat,S);
	return xopt;
}

void GWRDlg::SaveDataToTable()
{
	wxLogMessage("");
	if (!table_int) return;
	int n_obs = table_int->GetNumberRows();
	int time=0;
	std::vector<bool> save_undefs(n_obs);
	std::vector<double> beta(n_obs);
	std::vector<SaveToTableEntry> data(indep_mat.cols());

	wxString pre = "_est";
	wxString field_name;
	for (int i=0;i<data.size();i++)
	{
		if (i==0)
			field_name="Intercept";
		else 
			field_name=m_independentlist->GetString(i-1) + pre;

		for (int j=0;j<n_obs;j++)
		{
			beta[j] = outbetamat(j,i);
			save_undefs[j] = false;
		}
		data[i].d_val = &beta;
		data[i].undefined = &save_undefs;
		data[i].label = "Beta";
		data[i].field_default = field_name;
		data[i].type = GdaConst::double_type;


		int col = table_int->FindColId(field_name);
		if ( col == wxNOT_FOUND) {
			int col_insert_pos = table_int->GetNumberCols();
			int time_steps = 1;
			int m_length_val = GdaConst::default_dbf_long_len;
			int m_decimals_val = 0;

			if (data[i].type == GdaConst::double_type) {
				m_length_val = GdaConst::default_dbf_double_len;
				m_decimals_val = GdaConst::default_dbf_double_decimals;
			} else if (data[i].type == GdaConst::long64_type) {
				m_length_val = GdaConst::default_dbf_long_len;
				m_decimals_val = 0;
			} else if (data[i].type == GdaConst::string_type) {
				m_length_val = GdaConst::default_dbf_string_len;
				m_decimals_val = 0;
			}

			col = table_int->InsertCol(data[i].type, field_name, col_insert_pos, time_steps, m_length_val, m_decimals_val);

		}
		if (col > 0) {
			if (data[i].d_val) {
				table_int->SetColData(col, time, *data[i].d_val);
			} else if (data[i].l_val) {
				table_int->SetColData(col, time, *data[i].l_val);
			}
			if (data[i].undefined) {
				table_int->SetColUndefined(col, time, *data[i].undefined);
			}
		}
		
	}
	if (project->FindTableGrid()) project->FindTableGrid()->Refresh();
}

void GWRDlg::PrintGWRresult()
{
	vector<wxString> indep_name;
	wxString field_name;
	wxString pre="_est";
	wxString slog;
	logReport =wxEmptyString;
	int cnt = 0;

	slog << "------------------------------------------------";
	slog << "------------------------------------------------\n"; cnt++;
	slog << "		Results of Geographically ";
	slog << "Weighted Regression\n";cnt++;
	slog << "------------------------------------------------";
	slog << "------------------------------------------------\n"; cnt++;
	if(WeightModel==1)
		slog << "Kernel function:	gaussian \n"; cnt++;
	if (WeightModel==2)
		slog << "Kernel function:	bi-square \n"; cnt++;
	slog << wxString::Format("Bandwidth:   %7f\n",bw); cnt++;
	slog << "------------------------------";
	slog << "Summary of GWR coefficient estimates------------------------------\n";cnt++; 
	slog << "			   Min.          1st Qu.        Median        3rd Qu.        Max.\n";cnt++; 
	//slog << "\n\n\n";

	for (int i=0;i<outbetamat.cols();i++)
	{
		std::vector<double> vec(m_obs);
		if (i==0)
			field_name = "Intercept";
		else
			field_name = m_independentlist->GetString(i-1)+pre;
		indep_name.push_back(field_name);
		for (int j=0;j<m_obs;j++)
		{
			vec[j]=outbetamat(j,i);
		}
		sort(vec.begin(),vec.end());

		BreaksPoint(vec);
		slog << GenUtils::PadTrim(field_name, 18);
		slog << wxString::Format("  %12.5g   %12.6g   %12.6g   %12.6g   %12.6g\n",
			vec[0],breaks[0],breaks[1],breaks[2],vec[vec.size()-1]); cnt++;

	}
	slog << "------------------------------------------------";
	slog << "------------------------------------------------\n"; cnt++;
	slog << "All the GWR coefficients have been written into the data frame, i.e.\n";cnt++;
	for (int i=0;i<indep_name.size();i++)
	{
		slog << indep_name[i] + " , ";
	}
	slog << "\n";cnt++;
	slog << "-------------------------------------";
	slog << "Diagnostic information-------------------------------------\n";cnt++; 
	slog << wxString::Format(" Number of data points:   %5d\n",m_obs); cnt++;
	slog << wxString::Format(" AICc:  %12.6f\n",diagnosticVal[6]); cnt++;
	slog << wxString::Format(" AIC:  %12.6f\n",diagnosticVal[5]); cnt++;
	slog << wxString::Format(" RSS:  %12.6f\n",diagnosticVal[8]); cnt++;
	slog << wxString::Format(" R-square:  %12.6f\n",diagnosticVal[3]); cnt++;
	slog << wxString::Format(" Adjust R-square:  %12.6f\n",diagnosticVal[4]); cnt++;
	slog << "------------------------------------------------";
	slog << "------------------------------------------------\n"; cnt++;

	slog << "\n\n"; cnt++; cnt++;
	logReport << slog;
}

void GWRDlg::DisplayRegression(wxString dump)
{
	wxLogMessage("");
	wxDateTime now = wxDateTime::Now();
	dump = ">>" + now.FormatDate() + " " + now.FormatTime() +"\n\n"+ dump;

	if (regReportDlg == 0) {
		regReportDlg = new RegressionReportDlg(this, dump);
		regReportDlg->Connect(wxEVT_DESTROY, wxWindowDestroyEventHandler(GWRDlg::OnReportClose), NULL, this);
	} else {
		regReportDlg->AddNewReport(dump);
	}

	//regReportDlg = new RegressionReportDlg(this, dump);
	//regReportDlg->Connect(wxEVT_DESTROY, wxWindowDestroyEventHandler(GWRDlg::OnReportClose), NULL, this);
	regReportDlg->Show(true);
	regReportDlg->m_textbox->SetSelection(0, 0);
}
void GWRDlg::OnReportClose(wxWindowDestroyEvent& event)
{
	wxLogMessage("Close ReportDlg. ");
	regReportDlg = 0;
}

double GWRDlg::Percentile(double x, const std::vector<double>& v)
{
	int N = v.size();
	double Nd = (double) N;
	double p_0 = (100.0/Nd) * (1.0-0.5);
	double p_Nm1 = (100.0/Nd) * (Nd-0.5);
	if (x <= p_0) return v[0];
	if (x >= p_Nm1) return v[N-1];

	for (int i=1; i<N; i++) {
		double p_i = (100.0/Nd) * ((((double) i)+1.0)-0.5);
		if (x == p_i) return v[i];
		if (x < p_i) {
			double p_im1 = (100.0/Nd) * ((((double) i))-0.5);
			return v[i-1] + Nd*((x-p_im1)/100.0)*(v[i]-v[i-1]);
		}
	}
	return v[N-1]; // execution should never get here
}

void GWRDlg::BreaksPoint(const std::vector<double>& var)
{
	int num_cats = 4;
	breaks.resize(3);
	for (int i=0, iend=breaks.size(); i<iend; i++) {
		breaks[i] =
			Percentile(((i+1.0)*100.0)/((double) num_cats), var);
	}
}