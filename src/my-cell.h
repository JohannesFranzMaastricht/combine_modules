// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef MY_CELL_H_
#define MY_CELL_H_

#include "core/sim_object/sim_object.h"
#include "neuroscience/neuron_soma.h"
#include "neuroscience/event/new_neurite_extension_event.h"


namespace bdm {

// for EventHandler
using NewNeuriteExtensionEvent = bdm::experimental::neuroscience::NewNeuriteExtensionEvent;


BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
  BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, cell_type_, can_divide_, cell_color_);

 public:
  MyCellExt() {}
  explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {}

  /// If MyCell divides, daughter 2 copies the data members from the mother
  template <typename TMother>
  MyCellExt(const CellDivisionEvent& event, TMother* mother)
      : Base(event, mother) {
    cell_type_[kIdx] = mother->cell_type_[mother->kIdx];
    can_divide_[kIdx] = mother->can_divide_[mother->kIdx];
    cell_color_[kIdx] = mother->cell_color_[mother->kIdx];
  }

  /// If a cell divides, daughter keeps the same state from its mother.
  template <typename TDaughter>
  void EventHandler(const CellDivisionEvent& event, TDaughter* daughter) {
    Base::EventHandler(event, daughter);
  }

  /// \brief EventHandler to modify the data members of this soma
  /// after a new neurite extension event.
  ///
  /// Performs the transition mother to daughter 1
  /// \param event contains parameters for new neurite extension
  /// \param neurite pointer to new neurite
  /// \see NewNeuriteExtensionEvent
  template <typename TNeurite>
  void EventHandler(const NewNeuriteExtensionEvent& event, TNeurite* neurite) {
    // forward to experimental::neuroscience::NeuronSoma::EventHandler
    Base::EventHandler(event, neurite);
  }

  // getter and setter for our new data member
  void SetCellType(int t) { cell_type_[kIdx] = t; }
  int GetCellType() const { return cell_type_[kIdx]; }

  void SetCanDivide(int d) { can_divide_[kIdx] = d; }
  bool GetCanDivide() const { return can_divide_[kIdx]; }

  void SetCellColor(int cell_color) { cell_color_[kIdx] = cell_color; }
  int GetCellColor() const { return cell_color_[kIdx]; }

 private:
  // declare new data member and define their type
  // private data can only be accessed by public function and not directly
  vec<int> cell_type_;
  vec<bool> can_divide_;
  vec<int> cell_color_;
};

}  // namespace bdm

#endif  // MY_CELL_H_