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
#ifndef COMBINE_MODULES_H_
#define COMBINE_MODULES_H_

#include "biodynamo.h"
#include "biology-modules.h"
#include "core/substance_initializers.h"
#include "my-neurite.h"
#include "my-cell.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  experimental::neuroscience::InitModule();
  Simulation simulation(argc, argv);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();

  auto construct_soma = [](const std::array<double, 3> position) {
    auto* cell = new MyCell(position);
    cell->SetDiameter(6);
    cell->SetAdherence(0.4);
    cell->SetMass(1.0);
    cell->SetCanDivide(true);
    // Don't use gAllEventIds in the following line.
    // Otherwise this module will be copied to neurites during neurite extension.
    // This wil cause a fatal error because the GrowDivide module is only allowed
    // for simulation objects that derive from Cell. Neurites do not.
    cell->AddBiologyModule(new GrowDivide(8, 400, {CellDivisionEvent::kEventId})); //
    return cell;
  };

  auto add_dendrites = [](MyCell* soma) {
    MyNeurite my_neurite;
    auto* dendrite_apical = soma->ExtendNewNeurite({0, 0, 1}, &my_neurite)->SimObject::As<MyNeurite>();
    dendrite_apical->AddBiologyModule(new ApicalElongationBM());
    dendrite_apical->SetCanBranch(true);
    auto* dendrite_basal1 = soma->ExtendNewNeurite({0, 0, -1}, &my_neurite)->SimObject::As<MyNeurite>();
    dendrite_basal1->AddBiologyModule(new BasalElongationBM());
  };

// 1. Create cell
  auto* cell = construct_soma({0, 0, 0});
  rm->push_back(cell);

  // Add substances
  ModelInitializer::DefineSubstance(kSubstanceApical, "substance_apical", 0, 0,
                                    100);
  ModelInitializer::DefineSubstance(kSubstanceBasal, "substance_basal", 0, 0,
                                    100);
  ModelInitializer::InitializeSubstance(kSubstanceApical, "substance_apical", GaussianBand(param->max_bound_, 200, Axis::kZAxis));
  ModelInitializer::InitializeSubstance(kSubstanceApical, "substance_basal", GaussianBand(param->min_bound_, 200, Axis::kZAxis));
  // TODO initialize substances

// 2. Let it divide
  simulation.GetScheduler()->Simulate(10); // causes a bug in the visualization

// 3. Remove GrowDivide
  const auto& bms = cell->GetAllBiologyModules();
  for (size_t i = 0; i < bms.size(); i++) {
    auto* to_be_removed = dynamic_cast<GrowDivide*>(bms[i]);
    if (to_be_removed != nullptr) {
      cell->RemoveBiologyModule(to_be_removed);
      break;
    }
  }

// 4. Add dendrites
  add_dendrites(cell);

// 5. Let dendrites grow (add substances first)
  simulation.GetScheduler()->Simulate(10);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // COMBINE_MODULES_H_
