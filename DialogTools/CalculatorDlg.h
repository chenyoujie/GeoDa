/**
 * GeoDa TM, Copyright (C) 2011-2015 by Luc Anselin - all rights reserved
 *
 * This file is part of GeoDa.
 * 
 * GeoDa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoDa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GEODA_CENTER_CALCULATOR_DLG_H__
#define __GEODA_CENTER_CALCULATOR_DLG_H__

#include <map>
#include <set>
#include <vector>
#include <wx/button.h>
#include <wx/dialog.h>
//#include <wx/html/htmlwin.h>
#include <wx/webview.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/treectrl.h>
#include "../FramesManagerObserver.h"
#include "../DataViewer/TableStateObserver.h"
#include "../VarCalc/GdaParser.h"
#include "../GdaConst.h"

class FramesManager;
class TableState;
class TableInterface;
class Project;

class CalculatorDlg: public wxDialog,
public TableStateObserver
{    
public:
    CalculatorDlg(Project* project, wxWindow* parent,
				   const wxString& title = "Calculator", 
				   const wxPoint& pos = wxDefaultPosition,
				   const wxSize& size = wxDefaultSize,
				   long style = wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER );
	virtual ~CalculatorDlg();

	void OnClose(wxCloseEvent& e);
	void OnExprUpdate(wxCommandEvent& e);
	void OnTarget(wxCommandEvent& e);
	wxString h_title(const wxString& title, int level);
	void OnFuncHelpSel(wxTreeEvent& e);
	void OnFuncHelpDClick(wxTreeEvent& e);
	void OnNew(wxCommandEvent& e);
	void OnAssign(wxCommandEvent& e);
	void OnSelect(wxCommandEvent& e);
	void AssignOrSelect(bool assign);
		
	void ConnectToProject(Project* project);
	void DisconnectFromProject();

	/** Implementation of TableStateObserver interface */
	virtual void update(TableState* o);
	virtual bool AllowTimelineChanges() { return true; }
	virtual bool AllowGroupModify(const wxString& grp_nm) { return true; }
	virtual bool AllowObservationAddDelete() { return true; }
	
private:
	/** If enable false, target choices are cleared and and disabled.
	 If enable true and if table_int true, target choices are updated. */
	void InitTargetChoice(bool enable);
	void InitFuncHelpTree();
	
	/** Fill quick_parser_table with preview values for all numeric
	    fields in current project table */
	void UpdateQuickParserTable();
	/** Fill full_parser_table with values according to active_ident_list.
	 This operation involves copying full vectors of data for full evaluation
	 and should only be called from OnAssign */
	void UpdateFullParserTable();
	
	/** Checks that current expression is valid and sets variable
	 expr_valid to true iff valid. */
	void ValidateExpression();
	bool IsTargSel();
	
	bool all_init;
	
	Project* project;
	TableState* table_state;
	TableInterface* table_int;
	
	std::vector<GdaTokenDetails> tokens;
	std::map<wxString, GdaFVSmtPtr> quick_parser_table;
	std::map<wxString, GdaFVSmtPtr> full_parser_table;
	std::vector<GdaTokenDetails> eval_toks;
	std::set<wxString> active_ident_set;
	bool expr_valid;
	GdaLexer lexer;
	
	wxTreeCtrl* func_help_tree; // ID_FUNC_HELP_TREE
	//wxHtmlWindow* html_win; // ID_HTML_WIN
	wxWebView* html_win; // ID_HTML_WIN

	wxTextCtrl* expr_t_ctrl; // ID_EXPRESSION
	wxStaticText* msg_s_txt; // ID_MESSAGE
	wxStaticText* target_lbl_s_txt; // ID_TARGET_LBL
	wxStaticText* select_lbl_s_txt; // ID_SELECT_LBL
	wxChoice* target_choice; // ID_TARGET
	wxButton* new_btn; // ID_NEW
	wxButton* assign_btn; // ID_ASSIGN
	wxButton* select_btn; // ID_SELECT
	
	wxTextAttr attr_num;
	wxTextAttr attr_func;
	wxTextAttr attr_unknown_func;
	wxTextAttr attr_ident;
	wxTextAttr attr_unknown_ident;
};

#endif
