# -*- coding: UTF-8 -*-
#
# FileName     : parse_cohort
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-26 17:43
# Last Modified: 2024-01-26 17:43
# Modified By  : EastsunW
# -------------
# Description  : 解析队列配置文件，并从测序结果中获得样本信息
# -------------

import yaml
import re
from pathlib import Path
from natsort import natsorted

class Cohort:
    def __init__(self, cohort_yaml, reference_yaml, dependency_yaml):
        self.__parse_dependencies(dependency_yaml)
        self.__parse_cohort(cohort_yaml)
        self.__parse_references(reference_yaml)

    def __parse_cohort(self, cohort_yaml):
        with open(cohort_yaml) as f:
            cohort_config = yaml.load(f, Loader=yaml.FullLoader)
        # 队列名
        self.id = cohort_config["cohort_id"]
        # 参考基因组
        self.ref_build = cohort_config["reference_build"]
        # 解析目标
        self.__parse_targets(cohort_config)
        # 解析样本信息
        self.__parse_samples(cohort_config)

    def __parse_references(self, reference_yaml):
        """
        根据参考基因组配置解析参考信息，根据目标的不同，会需要不同的文件
        """
        with open(reference_yaml) as f:
            ref_config = yaml.load(f, Loader=yaml.FullLoader)
        if self.ref_build not in ref_config:
            raise KeyError(f"没有找到参考基因组{self.ref_build}")
        self.ref_config = ref_config[self.ref_build]

    def __parse_samples(self, config):
        """
        从配置文件的测序结果目录中，找到所有的不重复样本
        """
        self.regex = re.compile(config["regex"])
        self.input_dir = config["input_dir"]
        self.file_suffix = config["file_suffix"]
        # 遍历测序文件夹，找到所有的样本
        self.__samples = {}
        for file in Path(self.input_dir).glob(f'*.{self.file_suffix}'):
            file_match = self.regex.search(str(file.name))
            assert file_match, f"文件名不符合规范：{file}"
            sample, movie = file_match.group("sample"), \
                file_match.group("movie")
            if sample in self.__samples.keys():
                if movie in self.__samples[sample]:
                    print(f"样本{sample}的{movie}已经存在，将会被覆盖!")
                self.__samples[sample][movie] = str(file.name)
            else:
                self.__samples[sample] = {movie: str(file.name)}

    def __parse_dependencies(self, dependency_yaml):
        with open(dependency_yaml) as f:
            dependency_config = yaml.load(f, Loader=yaml.FullLoader)
        self.dependencies = dependency_config

    def __parse_targets(self, config):
        def __checker_recursive(requirements, steps_list, steps, missing_dependencies=None):
            if missing_dependencies is None:
                missing_dependencies = set()
            for step in steps_list:
                if step not in requirements or requirements.get(step) is None:
                    continue
                dependency = requirements[step]
                if dependency not in steps:
                    missing_dependencies.add(dependency)
                missing_dependencies = __checker_recursive(
                    requirements, [dependency], steps, missing_dependencies)
            return missing_dependencies
        requirements = self.dependencies
        steps = [step for step in config["targets"] if step in requirements]
        missing_dependencies = __checker_recursive(
            requirements, steps, steps)
        if missing_dependencies:
            print(f"自动添加缺失的依赖: {list(missing_dependencies)}")
            self.targets = steps + list(missing_dependencies)
        else:
            self.targets = steps

    def list_samples(self):
        return natsorted(list(self.__samples.keys()))

    def sample_info(self, sample):
        return self.__samples[sample]

    def list_movies(self, sample=None):
        if sample:
            return self.__samples[sample].keys()
        else:
            return [movie for sample in self.__samples for movie in self.__samples[sample].keys()]

    def list_files(self, sample=None, movie=None):
        if sample:
            if movie:
                return self.__samples[sample][movie]
            else:
                return [self.__samples[sample][movie] for movie in self.__samples[sample]]
        else:
            return [self.__samples[sample][movie] for sample in self.__samples for movie in self.__samples[sample]]

    def __str__(self):
        return (
            f'WGS Pacbio HiFi pipeline for HZAU Weibin Project\n' +
            f'cohort: {self.cohort}\n' +
            f'reference: {self.ref_build}\n' +
            f'targets: {self.targets}\n' +
            f'samples: {len(self.list_samples())}\n' +
            f'files: {len(self.list_files())}\n' +
            "\n".join(["  " + f'{sample}:\t{self.list_files(sample)}' for sample in self.list_samples()])
        )

    def __getitem__(self, item):
        if isinstance(item, str):
            return self.sample_info(item)
        if isinstance(item, int):
            name = self.list_samples()[item]
            return self.sample_info(name)
        raise TypeError("索引只能是整数或者字符串")

    def __len__(self):
        return len(self.list_samples())

    def __iter__(self):
        for sample in self.list_samples():
            yield sample

    def __contains__(self, item):
        return item in self.list_samples()

    def __repr__(self):
        return self.__str__()
