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

class MyCell : public experimental::neuroscience::NeuronSoma {
  BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, cell_type_, can_divide_, cell_color_);

 public:
  MyCell() : Base() {}

  virtual ~MyCell() {}

  explicit MyCell(const std::array<double, 3>& position) : Base(position) {}

  /// If MyCell divides, daughter 2 copies the data members from the mother
  MyCell(const Event& event, SimObject* other, uint64_t new_oid = 0)
      : Base(event, other, new_oid) {
    if (event.GetId() == CellDivisionEvent::kEventId) {
      auto* mother = static_cast<MyCell*>(other);
      cell_type_ = mother->cell_type_;
      can_divide_ = mother->can_divide_;
      cell_color_ = mother->cell_color_;
    }
  }

  void EventHandler(const Event& event, SimObject* other1,
                    SimObject* other2 = nullptr) {
    Base::EventHandler(event, other1, other2);
  }

  // getter and setter for our new data member
  void SetCellType(int t) { cell_type_ = t; }
  int GetCellType() const { return cell_type_; }

  void SetCanDivide(int d) { can_divide_ = d; }
  bool GetCanDivide() const { return can_divide_; }

  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  int GetCellColor() const { return cell_color_; }

 private:
  // declare new data member and define their type
  // private data can only be accessed by public function and not directly
  int cell_type_;
  bool can_divide_;
  int cell_color_;
};

}  // namespace bdm

#endif  // MY_CELL_H_
