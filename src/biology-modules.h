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
#ifndef BIOLOGY_MODULES_H_
#define BIOLOGY_MODULES_H_

#include "core/biology_module/biology_module.h"
#include "core/resource_manager.h"
#include "my-neurite.h"

namespace bdm {

enum Substances { kSubstanceApical, kSubstanceBasal };

struct ApicalElongationBM : public BaseBiologyModule {
  ApicalElongationBM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  ApicalElongationBM(const Event& event, BaseBiologyModule* other, uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new ApicalElongationBM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override { return new ApicalElongationBM(*this); }

  /// Default event handler (exising biology module won't be modified on
  /// any event)
  void EventHandler(const Event& event, BaseBiologyModule* other1,
                    BaseBiologyModule* other2 = nullptr) override {
    BaseBiologyModule::EventHandler(event, other1, other2);
  }

  // TODO: don't copy BM when split (elongate)

  void Run(SimObject* so) {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid(kSubstanceApical);
      init_ = true;
    }

    auto* dendrite = so->As<MyNeurite>();
    if (!dendrite) {
      return;
    }

    if (dendrite->GetDiameter() > 0.5) {
      std::array<double, 3> gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

      double gradient_weight = 0.04;
      double randomness_weight = 0.3;
      double old_direction_weight = 3;

      auto random_axis = random->template UniformArray<3>(-1, 1);

      auto old_direction =
          Math::ScalarMult(old_direction_weight, dendrite->GetSpringAxis());
      auto grad_direction = Math::ScalarMult(gradient_weight, gradient);
      auto random_direction = Math::ScalarMult(randomness_weight, random_axis);

      auto new_step_direction =
          Math::Add(Math::Add(old_direction, random_direction), grad_direction);

      dendrite->ElongateTerminalEnd(25, new_step_direction);
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.001);

      if (dendrite->GetCanBranch() && dendrite->IsTerminal() &&
          random->Uniform() < 0.033) {
        auto rand_noise = random->template UniformArray<3>(-0.1, 0.1);
        auto branch_direction = Math::Add(
            Math::Perp3(Math::Add(dendrite->GetUnitaryAxisDirectionVector(),
                                  rand_noise),
                        random->Uniform(0, 1)),
            dendrite->GetSpringAxis());
        auto* dendrite_2 = dendrite->Branch(branch_direction)->SimObject::As<MyNeurite>();
        dendrite_2->SetCanBranch(false);
        dendrite_2->SetDiameter(0.65);
      }

    }  // end if diameter
  }    // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  BDM_CLASS_DEF_OVERRIDE(ApicalElongationBM, 1);
};

struct BasalElongationBM : public BaseBiologyModule {
  BasalElongationBM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  BasalElongationBM(const Event& event, BaseBiologyModule* other,
                     uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new BasalElongationBM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new BasalElongationBM(*this);
  }

  void Run(SimObject* so) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid(kSubstanceBasal);
      init_ = true;
    }
    auto* dendrite = so->As<MyNeurite>();
    if (!dendrite) {
      return;
    }

    if (dendrite->IsTerminal() && dendrite->GetDiameter() > 0.75) {
      std::array<double, 3> gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

      double gradient_weight = 0.02;
      double randomness_weight = 0.5;
      double old_direction_weight = 5;

      auto random_axis = random->template UniformArray<3>(-1, 1);

      auto old_direction =
          Math::ScalarMult(old_direction_weight, dendrite->GetSpringAxis());
      auto grad_direction = Math::ScalarMult(gradient_weight, gradient);
      auto random_direction = Math::ScalarMult(randomness_weight, random_axis);

      auto new_step_direction =
          Math::Add(Math::Add(old_direction, random_direction), grad_direction);

      dendrite->ElongateTerminalEnd(25, new_step_direction);
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.001);

      if (random->Uniform() < 0.008) {
        dendrite->SetDiameter(dendrite->GetDiameter() - 0.01);
        dendrite->Bifurcate();
      }
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  BDM_CLASS_DEF_OVERRIDE(BasalElongationBM, 1);
};

}  // namespace bdm

#endif  // BIOLOGY_MODULES_H_
